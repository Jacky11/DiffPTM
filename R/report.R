#' Materials & Methods paragraph generated from the processing log
#'
#' Every function that changes the data appends one sentence to `x$log`; this function
#' assembles them with package versions so the text always matches what was run.
#'
#' @param x a `ptm_experiment`.
#' @param methods methods used (for citations).
#' @param blocking the `blocking` option of the usage analysis (`"auto"`, the per-site paired
#'   model, by default; `"random"` for the consensus single model; `"none"` for the unpaired one).
#' @param references append the full references cited in the paragraph.
#' @export
methods_paragraph <- function(x, methods = "limma", blocking = "auto", references = TRUE) {
  ver <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "n/a")
  cites <- c(limma = sprintf("limma v%s (Ritchie et al., 2015)", ver("limma")),
             msstatsptm = sprintf("MSstatsPTM v%s (Kohler et al., 2023)", ver("MSstatsPTM")),
             msqrob2ptm = sprintf("msqrob2PTM v%s (Demeulemeester et al., 2024)", ver("msqrob2")))
  imputed <- any(grepl("imp4p", x$log)); cp <- requireNamespace("cp4p", quietly = TRUE)
  txt <- paste0(
    input_files_sentence(x),
    "Data were processed with DiffPTM v", ver("DiffPTM"), " (R ", R.version$major, ".", R.version$minor,
    "), which extends the workflow of Giai Gianetto (2023) for the quality control, the detection profiles and the ",
    "differential analysis of PTMs relative to their protein across several biological conditions. ",
    paste(x$log, collapse = " "),
    design_sentence(x),
    " Differential analyses were performed with ", paste(cites[methods], collapse = ", "),
    ".",
    if ("limma" %in% methods) paste0(" In the single model (Giai Gianetto, 2023), site and matched protein intensities were fitted together (component x condition design)",
      if (!is.null(x$ptm$protein_intensity) && setequal(colnames(x$ptm$intensity), colnames(x$ptm$protein_intensity))) switch(blocking,
        random = "; the sample was included as a random effect with a consensus within-sample correlation (limma::duplicateCorrelation; consensus single model), and only samples in which both the site and its protein were quantified entered the usage model (complete pairs)",
        none = "; without a sample term (unpaired single model)",
        "; the sample was included as a fixed effect (per-site paired single model), so that each site's usage test is the paired analysis of its within-sample site-minus-protein differences") else "",
      "; differential usage was tested as the component x condition interaction with moderated t-statistics (robust empirical Bayes).") else "",
    if (imputed) paste0(" Missing values were imputed with imp4p v", ver("imp4p"), " (Giai Gianetto et al., 2020).") else "",
    if (cp) paste0(" P-values were adjusted with a pi0-aware FDR (cp4p v", ver("cp4p"), ", Giai Gianetto et al., 2016).") else " P-values were adjusted with the Benjamini-Hochberg procedure.")
  if (!references) return(txt)
  refs <- c(
    paste0("Giai Gianetto Q, Salhi Z. DiffPTM: differential abundance and usage analysis of post-translational modifications. ",
           "R package version ", ver("DiffPTM"), ". https://github.com/Jacky11/DiffPTM"),
    paste0("Giai Gianetto Q. Statistical Analysis of Post-Translational Modifications Quantified by Label-Free Proteomics Across ",
           "Multiple Biological Conditions with R: Illustration from SARS-CoV-2 Infected Cells. In: Burger T (ed.), Statistical ",
           "Analysis of Proteomic Data. Methods Mol Biol. 2023;2426:267-302. doi:10.1007/978-1-0716-1967-4_12"),
    if ("limma" %in% methods) "Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43:e47. doi:10.1093/nar/gkv007",
    if ("msstatsptm" %in% methods) "Kohler D, Tsai T-H, Verschueren E, et al. MSstatsPTM: statistical relative quantification of posttranslational modifications in bottom-up mass spectrometry-based proteomics. Mol Cell Proteomics. 2023;22:100477. doi:10.1016/j.mcpro.2022.100477",
    if ("msqrob2ptm" %in% methods) "Demeulemeester N, G\u00e9belin M, Caldi Gomes L, et al. msqrob2PTM: differential abundance and differential usage analysis of MS-based proteomics data at the posttranslational modification and peptidoform level. Mol Cell Proteomics. 2024;23:100708. doi:10.1016/j.mcpro.2023.100708",
    if (imputed) "Giai Gianetto Q, Wieczorek S, Cout\u00e9 Y, Burger T. A peptide-level multiple imputation strategy accounting for the different natures of missing values in proteomics data. bioRxiv. 2020. doi:10.1101/2020.05.29.122770",
    if (cp) "Giai Gianetto Q, Combes F, Ramus C, Bruley C, Cout\u00e9 Y, Burger T. Calibration plot for proteomics: a graphical tool to visually check the assumptions underlying FDR control in quantitative experiments. Proteomics. 2016;16:29-32. doi:10.1002/pmic.201500189")
  paste0(txt, "\n\nReferences\n", paste0("- ", refs, collapse = "\n"))
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
#' @param methods_text Materials & Methods text written to the README sheet (default:
#'   [methods_paragraph()]; the interface passes the text as edited by the user).
#' @details Every column of every sheet is explained in the `column_descriptions` sheet, and in
#'   a comment on its header cell in the result sheets.
#' @export
write_results_xlsx <- function(x, res, path, extra = list(), blocking = "auto", methods_text = NULL) {
  wb <- openxlsx::createWorkbook()
  desc <- list(); coldesc <- list()
  add <- function(name, df, description, samples = NULL, comments = FALSE) {
    name <- substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", name), 1, 31)
    openxlsx::addWorksheet(wb, name)
    df <- as.data.frame(df)
    openxlsx::writeDataTable(wb, name, df)
    desc[[length(desc) + 1]] <<- data.frame(sheet = name, description = description)
    cd <- describe_columns(names(df), samples)
    cdd <- data.frame(sheet = name, column = names(df), description = cd)
    if (!is.null(samples)) {                       # one row for all the sample columns
      isS <- cdd$column %in% names(samples)
      cdd <- rbind(cdd[!isS, ], data.frame(sheet = name, column = sprintf("one column per sample (%s%s)", paste(utils::head(names(samples), 3), collapse = ", "), if (length(samples) > 3) ", ..." else ""),
                                           description = sub("sample \\S+", "this sample", samples[[1]])))
    }
    coldesc[[length(coldesc) + 1]] <<- cdd
    if (comments) for (j in seq_along(cd)) if (nzchar(cd[j]))
      openxlsx::writeComment(wb, name, col = j, row = 1, comment = openxlsx::createComment(cd[j], author = "DiffPTM", width = 3, height = 4))
  }
  openxlsx::addWorksheet(wb, "README")
  for (a in c("protein", "ptm")) {
    if (is.null(x[[a]])) next
    smp <- colnames(x[[a]]$intensity)
    add(paste0(a, "_intensities"), cbind(x[[a]]$meta, x[[a]]$intensity),
        sprintf("Processed log2 intensities of the %s assay with feature metadata", a),
        samples = stats::setNames(sprintf("Processed log2 intensity in sample %s (after filtering and normalisation; empty = not observed)", smp), smp))
    add(paste0(a, "_detection"), cbind(feature_id = x[[a]]$meta$feature_id, qc_detection(x, a)),
        sprintf("Detection profile (1 = observed before imputation) of the %s assay", a),
        samples = stats::setNames(sprintf("1 if the feature was observed in sample %s (before any imputation), 0 if missing", smp), smp))
    add(paste0(a, "_jaccard"), cbind(sample = smp, as.data.frame(qc_jaccard(x, a))), "Jaccard index between samples on detection patterns",
        samples = stats::setNames(sprintf("Jaccard index with sample %s: features observed in both / features observed in either", smp), smp))
    add(paste0(a, "_pearson"), cbind(sample = smp, as.data.frame(qc_pearson(x, a))), "Pearson correlation between samples",
        samples = stats::setNames(sprintf("Pearson correlation of the log2 intensities with sample %s (features observed in both)", smp), smp))
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
    a <- if (an == "protein_DPA") "protein" else "ptm"                 # mean log2 intensity per group, as in the original app
    gm <- group_means(x[[a]]$intensity, groups_of(x, a))
    wide <- cbind(wide, gm[match(wide$feature_id, rownames(x[[a]]$intensity)), , drop = FALSE])
    if (an == "ptm_DPU" && !is.null(x$ptm$protein_intensity)) {
      pm <- group_means(x$ptm$protein_intensity, groups_of(x, "protein")[match(colnames(x$ptm$protein_intensity), x$annotation$protein$sample)])
      names(pm) <- paste0("protein_", names(pm))
      wide <- cbind(wide, pm[match(wide$feature_id, rownames(x$ptm$intensity)), , drop = FALSE])
    }
    what <- c(protein_DPA = "differential protein abundance", ptm_DPA = "differential PTM abundance (site intensity, not corrected for its protein)",
              ptm_DPU = "differential PTM usage (site corrected for its protein)")[[an]]
    add(paste(an, e, sep = "_"), wide, sprintf("Results of %s with method %s, one block of columns per contrast: log2FC, standard error, degrees of freedom, 95 %% confidence interval, p-value, adjusted p-value (NA = not tested). See column_descriptions.", what, e),
        comments = TRUE)
  }
  known <- c(detection_discordance = "PTM features detected in samples where their protein was not (and the reverse), on the detection pattern before imputation",
             method_comparison = "Significant calls and log2FC of every method side by side, per feature, analysis and contrast")
  for (nm in names(extra)) add(nm, extra[[nm]], if (nm %in% names(known)) known[[nm]] else paste("Additional table:", nm), comments = TRUE)
  cdf <- do.call(rbind, coldesc)
  openxlsx::addWorksheet(wb, "column_descriptions")
  openxlsx::writeDataTable(wb, "column_descriptions", cdf)
  openxlsx::setColWidths(wb, "column_descriptions", cols = 1:3, widths = c(28, 34, 120))
  if (is.null(methods_text)) methods_text <- methods_paragraph(x, unique(res$method), blocking)
  readme <- rbind(data.frame(sheet = "README", description = "Sheet descriptions and Materials & Methods"),
                  do.call(rbind, desc),
                  data.frame(sheet = "column_descriptions", description = "Meaning of every column of every sheet (also shown as a comment on the column headers of the result sheets)"),
                  data.frame(sheet = "Methods", description = methods_text))
  openxlsx::writeData(wb, "README", readme)
  openxlsx::setColWidths(wb, "README", cols = 1:2, widths = c(30, 120))
  openxlsx::addStyle(wb, "README", openxlsx::createStyle(wrapText = TRUE, valign = "top"), rows = seq_len(nrow(readme) + 1), cols = 2, gridExpand = TRUE)
  openxlsx::worksheetOrder(wb) <- c(1, length(openxlsx::sheets(wb)), seq(2, length(openxlsx::sheets(wb)) - 1))
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

#' Mean log2 intensity per group (observed values)
#' @keywords internal
group_means <- function(M, g) {
  g <- factor(g, levels = unique(as.character(g)))
  out <- vapply(levels(g), function(l) rowMeans(M[, g == l, drop = FALSE], na.rm = TRUE), numeric(nrow(M)))
  out <- matrix(out, nrow(M)); out[is.nan(out)] <- NA
  stats::setNames(as.data.frame(out), paste0("mean_log2_", levels(g)))
}

#' Plain-language description of result and table columns
#' @param cols column names.
#' @param samples optional named character vector: descriptions of sample columns.
#' @return character vector of descriptions (`""` when unknown).
#' @keywords internal
describe_columns <- function(cols, samples = NULL) {
  rules <- list(
    "^feature_id$" = "Identifier of the feature (protein group, PTM site or peptidoform); the same in every sheet",
    "^sample$" = "Sample",
    "^label$" = "Readable name of the feature (gene and site, e.g. GENE_S123, or protein)",
    "^gene$" = "Gene name(s), from the input table",
    "^protein_ids$" = "Protein accession(s), from the input table (several separated by ;)",
    "^stripped_sequence$" = "Peptide sequence without modifications",
    "^modified_sequence$" = "Peptide sequence with its modifications",
    "^residue$" = "Modified amino acid",
    "^position$" = "Position of the modified residue in the first protein",
    "^positions_within_proteins$" = "Position of the site in each protein of protein_ids",
    "^ptm_type$" = "Type of modification",
    "^localization_prob$" = "Localization probability of the site, from the software (1 = certain)",
    "^multiplicity$" = "Number of modifications on the measured peptide (MaxQuant ___1, ___2, ___3)",
    "^matched_protein$" = "feature_id of the proteome protein group(s) matched to this PTM feature; their abundance is used for the usage (DPU) correction",
    "^n_matched_protein$" = "Number of proteome protein groups matched to this PTM feature (0 = no usage test possible)",
    "^n_peptidoforms$" = "Number of peptidoforms summed into this site",
    "^multisite_key$" = "Key grouping the peptidoforms that carry the same set of sites",
    "^(Reverse|Potential\\.contaminant|Only\\.identified\\.by\\.site)$" = "Flag of the input table (TRUE = flagged by the software)",
    "^analysis$" = "protein_DPA = differential protein abundance; ptm_DPA = differential PTM abundance (site intensity); ptm_DPU = differential PTM usage (site corrected for its protein)",
    "^contrast$" = "Comparison: first group minus second group (B_vs_A = B minus A)",
    "^group$" = "Group (condition, and second factor if any)",
    "^n_samples$" = "Number of samples of the group",
    "^n_ptm_only$" = "Samples of the group where the PTM feature is observed and its protein is not",
    "^n_protein_only$" = "Samples of the group where the protein is observed and the PTM feature is not",
    "^n_both$" = "Samples of the group where both are observed",
    "^flag$" = "concordant; ptm_without_protein (the modification is seen where its protein is not: no usage test possible); protein_without_ptm; no_parent_protein (no matching protein in the proteome table)",
    "^n_methods_significant$" = "Number of methods calling the feature significant",
    "^protein_mean_log2_" = "Mean log2 intensity, in this group, of the protein matched to the PTM feature (observed values)",
    "^mean_log2_" = "Mean log2 intensity of the feature in this group (processed values, observed samples only)")
  out <- character(length(cols))
  for (i in seq_along(cols)) {
    v <- cols[i]
    if (!is.null(samples) && v %in% names(samples)) { out[i] <- samples[[v]]; next }
    hit <- which(vapply(names(rules), function(r) grepl(r, v, perl = TRUE), logical(1)))
    if (length(hit)) { out[i] <- rules[[hit[1]]]; next }
    m <- regmatches(v, regexec("^(log2FC|se|df|ci_lo|ci_hi|statistic|pvalue|adj_pvalue)[._](.+)$", v))[[1]]
    if (length(m)) {
      per <- if (substr(v, nchar(m[2]) + 1, nchar(m[2]) + 1) == ".") sprintf("contrast %s", m[3]) else sprintf("method %s", m[3])
      out[i] <- switch(m[2],
        log2FC = sprintf("log2 fold change, %s: first minus second group of the contrast; for usage (DPU), change of the site intensity relative to its protein (site minus protein)", per),
        se = sprintf("Standard error of the log2 fold change, %s", per),
        df = sprintf("Degrees of freedom of the test, %s (limma: residual plus prior degrees of freedom)", per),
        ci_lo = sprintf("Lower bound of the 95 %% confidence interval of the log2 fold change, %s", per),
        ci_hi = sprintf("Upper bound of the 95 %% confidence interval of the log2 fold change, %s", per),
        statistic = sprintf("Test statistic (moderated t), %s", per),
        pvalue = sprintf("Raw p-value, %s", per),
        adj_pvalue = sprintf("Adjusted p-value (false discovery rate), %s: the feature is significant when it is below the chosen FDR (e.g. 0.05)", per))
      next
    }
    m <- regmatches(v, regexec("^sig_(.+)$", v))[[1]]
    if (length(m)) { out[i] <- sprintf("TRUE if significant with method %s at the chosen FDR and fold-change thresholds", m[2]); next }
  }
  out
}

#' One sentence describing the groups of the design, for the Methods paragraph
#' @keywords internal
design_sentence <- function(x) {
  ann <- x$annotation$ptm %||% x$annotation$protein
  if (is.null(ann)) return("")
  cc <- unique(ann$condition)
  if ("factor" %in% names(ann) && !all(is.na(ann$factor)))
    sprintf(" Samples were annotated with a condition (%s) and a second factor (%s); contrasts compared conditions within each level of the second factor.",
            paste(cc, collapse = ", "), paste(unique(stats::na.omit(ann$factor)), collapse = ", "))
  else sprintf(" Samples were annotated with %d conditions (%s).", length(cc), paste(cc, collapse = ", "))
}

#' Sentence naming the input files, their size and modification date, for the Methods paragraph
#' @keywords internal
input_files_sentence <- function(x) {
  one <- function(src, what) {
    if (is.null(src) || !nrow(src)) return(NULL)
    sz <- src$size[1]
    size <- if (is.na(sz)) "" else if (sz >= 1024^2) sprintf("%.1f MB", sz / 1024^2) else sprintf("%.0f kB", sz / 1024)
    date <- if (is.na(src$modified[1])) "" else format(as.POSIXct(src$modified[1]), "%Y-%m-%d %H:%M")
    det <- paste(c(size, if (nzchar(date)) paste("modified", date)), collapse = ", ")
    sprintf("%s %s%s", what, src$name[1], if (nzchar(det)) sprintf(" (%s)", det) else "")
  }
  parts <- c(one(x$protein$source, "total proteome table"), one(x$ptm$source, "PTM table"), one(x$fasta_source, "FASTA"))
  if (!length(parts)) return("")
  src <- x$ptm$source %||% x$protein$source
  sw <- if (!is.null(src) && !is.na(src$software[1])) sprintf(", exported from %s", src$software[1]) else ""
  imp <- if (!is.null(src) && !is.null(src$imported)) sprintf(" and imported on %s", format(as.POSIXct(src$imported[1]), "%Y-%m-%d")) else ""
  sprintf("Input files%s%s: %s. ", sw, imp, paste(parts, collapse = "; "))
}
