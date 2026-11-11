mod_download_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("Download"),
    shiny::p(class = "step-help", "The Excel workbook contains all tables, a README sheet with the sheet descriptions and the Materials & Methods text below, ",
             "and a column_descriptions sheet explaining every column (also shown as comments on the column headers). ",
             "Figures are exported as editable PowerPoint and as PDF: quality control (correlation and detection heatmaps, PCA, ",
             "intensity distributions) first, then the results. Building the figure files can take a few seconds."),
    shiny::downloadButton(ns("xlsx"), "Excel workbook"),
    shiny::downloadButton(ns("pptx"), "PowerPoint (editable figures)"),
    shiny::downloadButton(ns("pdf"), "PDF figures"),
    shiny::downloadButton(ns("rds"), "Session (.rds, reload in R with readRDS)"),
    shiny::h4("Materials & Methods"),
    shiny::p(class = "step-help", "Written from what was run. You can edit it here before downloading the Excel workbook (the edited text is the one exported)."),
    shiny::uiOutput(ns("mm_ui")),
    shiny::actionButton(ns("mm_reset"), "Restore the generated text", icon = shiny::icon("rotate-left"))
  )
}

mod_download_server <- function(id, stats) {
  shiny::moduleServer(id, function(input, output, session) {
    generated <- shiny::reactive({ s <- stats(); shiny::req(s); methods_paragraph(s$x, s$methods, if (is.null(s$blocking)) "auto" else s$blocking) })
    reset <- shiny::reactiveVal(0)
    shiny::observeEvent(input$mm_reset, reset(reset() + 1))
    output$mm_ui <- shiny::renderUI({ reset()
      shiny::textAreaInput(session$ns("mm_text"), NULL, generated(), width = "100%", rows = 12, resize = "vertical") })
    shiny::outputOptions(output, "mm_ui", suspendWhenHidden = FALSE)
    mm <- function() { t <- input$mm_text; if (is.null(t) || !nzchar(trimws(t))) generated() else t }
    figures <- function(s) {
      figs <- list()
      nm <- c(protein = "total proteome", ptm = "PTM")
      safe <- function(e) tryCatch(e, error = function(err) NULL)
      # 1. heatmaps and clustering, 2. detection overlaps, 3. other quality control, 4. results
      for (a in c("protein", "ptm")) {
        figs[[paste("Correlation matrix,", nm[[a]])]] <- safe(plot_correlation_matrix(s$x, a))
        figs[[paste("Jaccard index on detection,", nm[[a]])]] <- safe(plot_jaccard_matrix(s$x, a))
      }
      for (a in c("protein", "ptm")) {
        figs[[paste("Clustering (Pearson),", nm[[a]])]] <- safe(plot_dendrogram(s$x, a, "pearson"))
        figs[[paste("Clustering (Jaccard),", nm[[a]])]] <- safe(plot_dendrogram(s$x, a, "jaccard"))
      }
      lv <- list(protein = list(c("protein", "feature", "proteins")), ptm = list(c("ptm", "feature", "PTM features"), c("ptm", "protein", "modified proteins")))
      for (L in unlist(lv, recursive = FALSE)) {
        figs[[paste("Detected", L[3], "- all samples")]] <- safe(plot_detection_overlap(s$x, L[1], level = L[2]))
        for (g in levels(groups_of(s$x, L[1]))) figs[[paste("Detected", L[3], "- condition", g)]] <- safe(plot_detection_overlap(s$x, L[1], group = g, level = L[2]))
      }
      for (a in c("protein", "ptm")) {
        figs[[paste("PCA,", nm[[a]])]] <- safe(plot_pca(s$x, a))
        figs[[paste("PCA variance explained,", nm[[a]])]] <- safe(plot_pca_scree(s$x, a))
        figs[[paste("Intensity distributions,", nm[[a]])]] <- safe(plot_intensity_distributions(s$x, a))
        figs[[paste("Features per sample,", nm[[a]])]] <- safe(plot_detected_per_sample(s$x, a))
      }
      lab <- function(r, an) { meta <- if (an == "protein_DPA") s$x$protein$meta else s$x$ptm$meta; r$label <- meta$label[match(r$feature_id, meta$feature_id)]; r }
      for (an in unique(s$res$analysis)) figs[[paste("Counts", an)]] <- plot_counts(s$res, s$fdr, s$min_log2fc, an)
      if (length(s$methods) > 1) for (an in unique(s$res$analysis)) for (ct in s$contrasts$name)
        figs[[paste("Agreement between methods", an, ct)]] <- plot_method_overlap(s$res, s$fdr, s$min_log2fc, an, ct)
      for (e in s$methods) for (ct in s$contrasts$name) {
        for (an in unique(s$res$analysis)) {
          r <- s$res[s$res$method == e & s$res$analysis == an & s$res$contrast == ct, ]
          if (!nrow(r)) next
          figs[[paste("Volcano", an, ct, e)]] <- plot_volcano(lab(r, an), s$fdr, s$min_log2fc, title = paste(an, ct, e))
          figs[[paste("p-values", an, ct, e)]] <- plot_pvalue_histogram(r, title = paste("Raw p-values", an, ct, e))
          cr <- s$contrasts[s$contrasts$name == ct, ]
          if (an != "ptm_DPU") figs[[paste("Mean intensities", an, ct, e)]] <- safe(plot_condition_means(s$x, r, cr$group1, cr$group2, if (an == "protein_DPA") "protein" else "ptm",
                                                                                     s$fdr, s$min_log2fc, title = paste("Mean intensities,", an, ct, e)))
        }
        r <- s$res[s$res$method == e & s$res$contrast == ct, ]
        if (all(c("ptm_DPA", "protein_DPA") %in% r$analysis)) {
          parent <- stats::setNames(s$x$ptm$meta$matched_protein, s$x$ptm$meta$feature_id)
          figs[[paste("PTM vs protein", ct, e)]] <- plot_fc_comparison(lab(r[r$analysis == "ptm_DPA", ], "ptm_DPA"), r[r$analysis == "protein_DPA", ], parent,
                                                                         r[r$analysis == "ptm_DPU", ], s$fdr, s$min_log2fc, title = paste(ct, e))
        }
      }
      figs[!vapply(figs, is.null, TRUE)]
    }
    with_progress <- function(msg, expr) shiny::withProgress(message = msg, value = 0.3, expr)
    output$xlsx <- shiny::downloadHandler("DiffPTM_results.xlsx", function(f) { s <- stats(); shiny::req(s)
      with_progress("Writing the Excel workbook", {
        extra <- if (length(s$methods) > 1) list(method_comparison = compare_methods(s$res, s$fdr, s$min_log2fc)$wide) else list()
        extra$detection_discordance <- detection_discordance(s$x)
        write_results_xlsx(s$x, s$res, f, extra, if (is.null(s$blocking)) "auto" else s$blocking, methods_text = mm()) }) })
    output$pptx <- shiny::downloadHandler("DiffPTM_figures.pptx", function(f) { s <- stats(); shiny::req(s)
      with_progress("Building the PowerPoint figures", write_results_pptx(figures(s), f)) })
    output$pdf  <- shiny::downloadHandler("DiffPTM_figures.pdf", function(f) { s <- stats(); shiny::req(s)
      with_progress("Building the PDF figures", write_results_pdf(figures(s), f)) })
    output$rds  <- shiny::downloadHandler("DiffPTM_session.rds", function(f) saveRDS(stats(), f))
    for (o in c("xlsx", "pptx", "pdf", "rds")) shiny::outputOptions(output, o, suspendWhenHidden = FALSE)
  })
}
