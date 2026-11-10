mod_results_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("Results"),
    shiny::fluidRow(
      shiny::column(2,
        shiny::selectInput(ns("analysis"), "Analysis", c("PTM usage, corrected for protein (DPU)" = "ptm_DPU", "PTM abundance (DPA)" = "ptm_DPA", "Protein abundance (DPA)" = "protein_DPA")),
        shiny::uiOutput(ns("contrast_ui")),
        shiny::uiOutput(ns("method_ui"))),
      shiny::column(5, shiny::plotOutput(ns("volcano"), height = 520)),
      shiny::column(5, shiny::plotOutput(ns("fc"), height = 520))
    ),
    shiny::fluidRow(
      shiny::column(4, shiny::plotOutput(ns("hist"), height = 300)),
      shiny::column(4, shiny::plotOutput(ns("counts"), height = 300)),
      shiny::column(4, shiny::h4("Agreement between methods"),
                    shiny::p(class = "step-help", "Significant features of the selected analysis and contrast, by combination of methods (UpSet plot: ",
                             "each bar counts the features called by exactly the methods marked below it). Run two or three methods at the Statistics step to see it."),
                    shiny::plotOutput(ns("overlap"), height = 330))
    ),
    shiny::h4("Significant features"), DT::DTOutput(ns("table")),
    shiny::h4("Detection discordance"),
    shiny::p(class = "step-help", "PTM features observed in at least two samples of a condition in which their protein was not quantified (ptm_without_protein; no usage test is possible, but the modification is there), and the reverse (protein_without_ptm). Computed on the detection pattern before any imputation (as in Giai Gianetto, Methods Mol Biol 2023, section 3.8)."),
    shiny::fluidRow(shiny::column(3, shiny::selectInput(ns("disc_flag"), "Pattern", c("PTM without protein" = "ptm_without_protein", "Protein without PTM" = "protein_without_ptm", "PTM without any protein in the table" = "no_parent_protein")))),
    DT::DTOutput(ns("discordance"))
  )
}

mod_results_server <- function(id, stats) {
  shiny::moduleServer(id, function(input, output, session) {
    output$contrast_ui <- shiny::renderUI({ s <- stats(); shiny::req(s); shiny::selectInput(session$ns("contrast"), "Contrast", s$contrasts$name) })
    output$method_ui   <- shiny::renderUI({ s <- stats(); shiny::req(s); shiny::selectInput(session$ns("method"), "Method", s$methods) })
    sub <- shiny::reactive({ s <- stats(); shiny::req(s, input$contrast, input$method)
      r <- s$res[s$res$analysis == input$analysis & s$res$contrast == input$contrast & s$res$method == input$method, ]
      meta <- if (input$analysis == "protein_DPA") s$x$protein$meta else s$x$ptm$meta
      cbind(r, meta[match(r$feature_id, meta$feature_id), setdiff(names(meta), "feature_id"), drop = FALSE]) })
    titles <- c(ptm_DPU = "PTM usage (DPU)", ptm_DPA = "PTM abundance (DPA)", protein_DPA = "Protein abundance (DPA)")
    output$volcano <- shiny::renderPlot({ s <- stats(); plot_volcano(sub(), s$fdr, s$min_log2fc, title = paste0(titles[input$analysis], ": ", input$contrast, " (", input$method, ")")) })
    output$hist    <- shiny::renderPlot(plot_pvalue_histogram(sub(), title = "Raw p-values"))
    output$fc <- shiny::renderPlot({ s <- stats(); shiny::req(input$analysis != "protein_DPA")
      r <- s$res[s$res$contrast == input$contrast & s$res$method == input$method, ]
      meta <- s$x$ptm$meta; r$label <- meta$label[match(r$feature_id, meta$feature_id)]
      parent <- stats::setNames(meta$matched_protein, meta$feature_id)
      plot_fc_comparison(r[r$analysis == "ptm_DPA", ], r[r$analysis == "protein_DPA", ], parent, r[r$analysis == "ptm_DPU", ], s$fdr, s$min_log2fc,
                         title = "PTM feature vs protein fold change") })
    output$counts <- shiny::renderPlot({ s <- stats(); plot_counts(s$res, s$fdr, s$min_log2fc, input$analysis) + ggplot2::labs(title = paste("Significant features,", titles[input$analysis])) })
    output$table <- DT::renderDT({ s <- stats(); d <- sub(); d <- d[!is.na(d$adj_pvalue) & d$adj_pvalue < s$fdr & abs(d$log2FC) >= s$min_log2fc, ]
      first <- intersect(c("label", "feature_id", "gene", "protein_ids", "residue", "position", "contrast", "log2FC", "se", "df", "ci_lo", "ci_hi", "statistic", "pvalue", "adj_pvalue"), names(d))
      d <- d[order(d$pvalue), c(first, setdiff(names(d), c(first, "analysis", "method")))]
      DT::datatable(d, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE) |> DT::formatSignif(intersect(c("log2FC", "se", "df", "ci_lo", "ci_hi", "statistic", "pvalue", "adj_pvalue"), names(d)), 3) })
    output$overlap <- shiny::renderPlot({ s <- stats(); shiny::req(length(s$methods) > 1, input$contrast)
      plot_method_overlap(s$res, s$fdr, s$min_log2fc, input$analysis, input$contrast, title = NULL, subtitle = NULL) +
        ggplot2::labs(title = NULL) })
    output$discordance <- DT::renderDT({ s <- stats(); d <- detection_discordance(s$x); d <- d[d$flag == input$disc_flag, ]
      meta <- s$x$ptm$meta; d <- cbind(label = meta$label[match(d$feature_id, meta$feature_id)], d)
      DT::datatable(d, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE) })
  })
}
