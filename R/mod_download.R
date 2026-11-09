mod_download_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("Download"),
    shiny::p(class = "step-help", "The Excel workbook contains all tables, a README sheet with sheet descriptions and the auto-generated Materials & Methods paragraph. Figures are exported as editable PowerPoint and PDF."),
    shiny::downloadButton(ns("xlsx"), "Excel workbook"),
    shiny::downloadButton(ns("pptx"), "PowerPoint (editable figures)"),
    shiny::downloadButton(ns("pdf"), "PDF figures"),
    shiny::downloadButton(ns("rds"), "Session (.rds, reload in R with readRDS)"),
    shiny::h4("Materials & Methods"), shiny::verbatimTextOutput(ns("mm_text"))
  )
}

mod_download_server <- function(id, stats) {
  shiny::moduleServer(id, function(input, output, session) {
    figures <- shiny::reactive({ s <- stats(); shiny::req(s)
      figs <- list()
      lab <- function(r, an) { meta <- if (an == "protein_DPA") s$x$protein$meta else s$x$ptm$meta; r$label <- meta$label[match(r$feature_id, meta$feature_id)]; r }
      for (an in unique(s$res$analysis)) figs[[paste("Counts", an)]] <- plot_counts(s$res, s$fdr, s$min_log2fc, an)
      for (e in s$methods) for (ct in s$contrasts$name) {
        for (an in unique(s$res$analysis)) {
          r <- s$res[s$res$method == e & s$res$analysis == an & s$res$contrast == ct, ]
          if (!nrow(r)) next
          figs[[paste("Volcano", an, ct, e)]] <- plot_volcano(lab(r, an), s$fdr, s$min_log2fc, title = paste(an, ct, e))
          figs[[paste("p-values", an, ct, e)]] <- plot_pvalue_histogram(r, title = paste("Raw p-values", an, ct, e))
        }
        r <- s$res[s$res$method == e & s$res$contrast == ct, ]
        if (all(c("ptm_DPA", "protein_DPA") %in% r$analysis)) {
          parent <- stats::setNames(s$x$ptm$meta$matched_protein, s$x$ptm$meta$feature_id)
          figs[[paste("PTM vs protein", ct, e)]] <- plot_fc_comparison(lab(r[r$analysis == "ptm_DPA", ], "ptm_DPA"), r[r$analysis == "protein_DPA", ], parent,
                                                                         r[r$analysis == "ptm_DPU", ], s$fdr, s$min_log2fc, title = paste(ct, e))
        }
      }
      for (a in c("protein", "ptm")) {
        figs[[paste("Detected", a)]] <- plot_detected_per_sample(s$x, a)
        figs[[paste("Jaccard", a)]] <- plot_heatmap_matrix(qc_jaccard(s$x, a), paste("Jaccard index on detection,", a),
                                                           groups = stats::setNames(as.character(groups_of(s$x, a)), s$x$annotation[[a]]$sample))
        figs[[paste("PCA", a)]] <- tryCatch(plot_pca(s$x, a), error = function(e) NULL)
        figs[[paste("Intensities", a)]] <- plot_intensity_distributions(s$x, a)
      }
      figs[!vapply(figs, is.null, TRUE)] })
    output$mm_text <- shiny::renderText({ s <- stats(); shiny::req(s); methods_paragraph(s$x, s$methods, if (is.null(s$blocking)) "auto" else s$blocking) })
    output$xlsx <- shiny::downloadHandler("DiffPTM_results.xlsx", function(f) { s <- stats()
      extra <- if (length(s$methods) > 1) list(method_comparison = compare_methods(s$res, s$fdr, s$min_log2fc)$wide) else list()
      extra$detection_discordance <- detection_discordance(s$x)
      write_results_xlsx(s$x, s$res, f, extra, if (is.null(s$blocking)) "auto" else s$blocking) })
    output$pptx <- shiny::downloadHandler("DiffPTM_figures.pptx", function(f) write_results_pptx(figures(), f))
    output$pdf  <- shiny::downloadHandler("DiffPTM_figures.pdf",  function(f) write_results_pdf(figures(), f))
    output$rds  <- shiny::downloadHandler("DiffPTM_session.rds",  function(f) saveRDS(stats(), f))
  })
}
