#' Materials & Methods paragraph generated from the processing log
#'
#' Every function that changes the data appends one sentence to `x$log`; this function
#' assembles them with package versions so the text always matches what was run.
#'
#' @param x a `ptm_experiment`.
#' @param methods methods used (for citations).
#' @param blocking the `blocking` option of the usage analysis (`"auto"`, the per-site paired
#'   model, by default; `"random"` for the consensus single model; `"none"` for the unpaired one).
#' @export
methods_paragraph <- function(x, methods = "limma", blocking = "auto") {
  ver <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "n/a")
  cites <- c(limma = sprintf("limma v%s (Ritchie et al., 2015)", ver("limma")),
             msstatsptm = sprintf("MSstatsPTM v%s (Kohler et al., 2023)", ver("MSstatsPTM")),
             msqrob2ptm = sprintf("msqrob2PTM v%s (Demeulemeester et al., 2024)", ver("msqrob2")))
  paste0(
    "Data were processed with DiffPTM v", ver("DiffPTM"), " (R ", R.version$major, ".", R.version$minor, "). ",
    paste(x$log, collapse = " "),
    " Differential analyses were performed with ", paste(cites[methods], collapse = ", "),
    ".",
    if ("limma" %in% methods) paste0(" In the single model, site and matched protein intensities were fitted together (component x condition design)",
      if (!is.null(x$ptm$protein_intensity) && setequal(colnames(x$ptm$intensity), colnames(x$ptm$protein_intensity))) switch(blocking,
        random = "; the sample was included as a random effect with a consensus within-sample correlation (limma::duplicateCorrelation; consensus single model), and only samples in which both the site and its protein were quantified entered the usage model (complete pairs)",
        none = "; without a sample term (unpaired single model)",
        "; the sample was included as a fixed effect (per-site paired single model), so that each site's usage test is the paired analysis of its within-sample site-minus-protein differences") else "",
      "; differential usage was tested as the component x condition interaction with moderated t-statistics (robust empirical Bayes).") else "",
    if (any(grepl("imp4p", x$log))) paste0(" Missing values were imputed with imp4p v", ver("imp4p"), " (Giai Gianetto et al., 2020).") else "",
    if (requireNamespace("cp4p", quietly = TRUE)) paste0(" P-values were adjusted with a pi0-aware FDR (cp4p v", ver("cp4p"), ", Giai Gianetto et al., 2016).") else " P-values were adjusted with the Benjamini-Hochberg procedure.")
}

#' Write results and QC tables to an Excel workbook
#'
#' One sheet per table plus a `README` sheet with sheet descriptions and the Methods
#' paragraph (replaces the `*_descriptions.txt` files of the original app).
#'
#' @param x a `ptm_experiment`.
#' @param res tidy results from [run_methods()].
#' @param path output `.xlsx`.
#' @param extra named list of additional data.frames (e.g. `compare_methods()$wide`).
#' @param blocking the pairing option of the usage analysis, for the Methods paragraph (see [methods_paragraph()]).
#' @export
write_results_xlsx <- function(x, res, path, extra = list(), blocking = "auto") {
  wb <- openxlsx::createWorkbook()
  desc <- list()
  add <- function(name, df, description) {
    name <- substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", name), 1, 31)
    openxlsx::addWorksheet(wb, name)
    openxlsx::writeDataTable(wb, name, as.data.frame(df))
    desc[[length(desc) + 1]] <<- data.frame(sheet = name, description = description)
  }
  openxlsx::addWorksheet(wb, "README")
  for (a in c("protein", "ptm")) {
    if (is.null(x[[a]])) next
    add(paste0(a, "_intensities"), cbind(x[[a]]$meta, x[[a]]$intensity),
        sprintf("Processed log2 intensities of the %s assay with feature metadata", a))
    add(paste0(a, "_detection"), cbind(feature_id = x[[a]]$meta$feature_id, qc_detection(x, a)),
        sprintf("Detection profile (1 = observed before imputation) of the %s assay", a))
    add(paste0(a, "_jaccard"), as.data.frame(qc_jaccard(x, a)), "Jaccard index between samples on detection patterns")
    add(paste0(a, "_pearson"), as.data.frame(qc_pearson(x, a)), "Pearson correlation between samples")
  }
  for (an in unique(res$analysis)) for (e in unique(res$method)) {
    sub <- res[res$analysis == an & res$method == e, ]
    if (!nrow(sub)) next
    cols <- intersect(c("feature_id", "contrast", "log2FC", "se", "df", "ci_lo", "ci_hi", "pvalue", "adj_pvalue"), names(sub))
    wide <- stats::reshape(sub[, cols], idvar = "feature_id", timevar = "contrast", direction = "wide")
    meta <- if (an == "protein_DPA") x$protein$meta else x$ptm$meta
    m <- meta[match(wide$feature_id, meta$feature_id), , drop = FALSE]
    first <- intersect(c("feature_id", "label", "gene", "protein_ids", "residue", "position", "multiplicity", "localization_prob"), names(m))
    wide <- cbind(m[, c(first, setdiff(names(m), c(first, "feature_id"))), drop = FALSE], wide[, setdiff(names(wide), "feature_id"), drop = FALSE])
    add(paste(an, e, sep = "_"), wide, sprintf("%s results of method %s per contrast: log2FC, standard error, degrees of freedom, 95 %% confidence interval, p-value, adjusted p-value (NA = not tested by this method)", an, e))
  }
  known <- c(detection_discordance = "PTM features detected in samples where their protein was not (and the reverse), on the detection pattern before imputation",
             method_comparison = "Significant calls and log2FC of every method side by side, per feature, analysis and contrast")
  for (nm in names(extra)) add(nm, extra[[nm]], if (nm %in% names(known)) known[[nm]] else paste("Additional table:", nm))
  readme <- rbind(data.frame(sheet = "README", description = "Sheet descriptions and Materials & Methods"),
                  do.call(rbind, desc),
                  data.frame(sheet = "Methods", description = methods_paragraph(x, unique(res$method), blocking)))
  openxlsx::writeData(wb, "README", readme)
  openxlsx::setColWidths(wb, "README", cols = 1:2, widths = c(30, 120))
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  invisible(path)
}

#' Write figures to an editable PowerPoint (one slide per figure, vector graphics)
#'
#' @param figures named list of ggplot objects (names become slide titles).
#' @param path output file.
#' @param title title slide text.
#' @export
write_results_pptx <- function(figures, path, title = "DiffPTM report") {
  doc <- officer::read_pptx()
  doc <- officer::add_slide(doc, layout = "Title Slide", master = "Office Theme")
  doc <- officer::ph_with(doc, title, location = officer::ph_location_type(type = "ctrTitle"))
  for (nm in names(figures)) {
    doc <- officer::add_slide(doc, layout = "Title and Content", master = "Office Theme")
    doc <- officer::ph_with(doc, nm, location = officer::ph_location_type(type = "title"))
    if (requireNamespace("rvg", quietly = TRUE)) {
      doc <- officer::ph_with(doc, rvg::dml(ggobj = figures[[nm]]), location = officer::ph_location_type(type = "body"))
    } else {                       # fallback: raster image when rvg is not installed
      png <- tempfile(fileext = ".png"); ggplot2::ggsave(png, figures[[nm]], width = 8, height = 6, dpi = 150)
      doc <- officer::ph_with(doc, officer::external_img(png), location = officer::ph_location_type(type = "body"))
    }
  }
  print(doc, target = path)
  invisible(path)
}

#' Write figures to a PDF (one page per figure)
#' @inheritParams write_results_pptx
#' @param width,height page size in inches.
#' @export
write_results_pdf <- function(figures, path, width = 8, height = 7) {
  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off())
  for (nm in names(figures)) print(figures[[nm]])
  invisible(path)
}
