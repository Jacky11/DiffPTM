mod_process_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("Filtering, matching, normalisation, imputation"),
    shiny::p(class = "step-help", "Every option chosen here is written to the Methods paragraph of the report."),
    shiny::fluidRow(
      shiny::column(4,
        shiny::h4("Filtering"),
        shiny::checkboxInput(ns("flags"), "Remove reverse hits / contaminants", TRUE),
        shiny::checkboxInput(ns("merge_redundant"), "Merge site rows repeated for several protein accessions (identical quantities)", TRUE),
        shiny::numericInput(ns("min_prob"), "Minimal localization probability", 0.75, 0, 1, 0.05),
        shiny::numericInput(ns("min_n"), "Min. observed values in at least one group", 2, 0),
        shiny::numericInput(ns("min_pct"), "Min. % observed values in at least one group", 0, 0, 100)
      ),
      shiny::column(4,
        shiny::h4("PTM to protein matching"),
        shiny::radioButtons(ns("match"), NULL, c("By protein identifier" = "id", "By peptide sequence (FASTA)" = "sequence", "Hybrid" = "hybrid")),
        shiny::fileInput(ns("fasta"), "FASTA of the search (sequence / hybrid matching; site positions of peptide tables)", accept = c(".fasta", ".fa", ".fas", ".txt")),
        shiny::radioButtons(ns("shared"), "Peptides shared by several proteins", c("Average parent proteins" = "mean", "First protein" = "first", "Discard" = "drop"), inline = TRUE),
        shiny::checkboxInput(ns("sites"), "Collapse peptidoforms to sites (sum of the peptidoforms of a site)", FALSE)
      ),
      shiny::column(4,
        shiny::h4("Normalisation"),
        shiny::radioButtons(ns("norm"), NULL, c("Global median centering (default)" = "global_median", "Per-condition median centering (keeps shifts between conditions)" = "median", "ComBat (second factor = batch)" = "combat", "None" = "none")),
        shiny::h4("Imputation"),
        shiny::radioButtons(ns("impute"), NULL, c("None (recommended for DPU)" = "none", "imp4p MLE" = "mle", "imp4p SLSA" = "slsa", "imp4p multiple imputation" = "mi")),
        shiny::actionButton(ns("run"), "Run processing", class = "btn-primary", icon = shiny::icon("gears"))
      )
    ),
    shiny::hr(),
    shiny::h4("Quality control"),
    shiny::fluidRow(
      shiny::column(3, shiny::selectInput(ns("qc_assay"), "Assay", c("PTM" = "ptm", "Total proteome" = "protein"))),
      shiny::column(9, shiny::verbatimTextOutput(ns("log")))
    ),
    shiny::fluidRow(
      shiny::column(6, shiny::plotOutput(ns("detected"), height = 300), shiny::plotOutput(ns("dist"), height = 300)),
      shiny::column(6, shiny::plotOutput(ns("jaccard"), height = 420), shiny::plotOutput(ns("pca"), height = 380))
    )
  )
}

mod_process_server <- function(id, annotated) {
  shiny::moduleServer(id, function(input, output, session) {
    x <- shiny::eventReactive(input$run, {
      x <- annotated(); shiny::req(x)
      shiny::withProgress(message = "Processing", {
        tryCatch({
          if (input$flags) x <- filter_flagged(x)
          if (input$merge_redundant) x <- tryCatch(merge_redundant_ptm(x), error = function(e) x)
          if (all(is.na(x$ptm$meta$position)) && !all(is.na(x$ptm$meta$modified_sequence)))     # peptidoform tables (DIA-NN, FragPipe modified peptides)
            x <- tryCatch(locate_sites(x, input$fasta$datapath), error = function(e) x)
          x <- filter_localization(x, input$min_prob)
          x <- filter_min_observed(x, min_n = input$min_n, min_pct = input$min_pct)
          shiny::incProgress(0.2, "Normalisation")
          x <- normalise(x, method = input$norm)
          shiny::incProgress(0.2, "Imputation")
          x <- impute(x, method = input$impute)
          shiny::incProgress(0.2, "Matching")
          if (input$sites) x <- collapse_to_sites(x)
          x <- match_ptm_to_protein(x, method = input$match, fasta = input$fasta$datapath, shared = input$shared)
          x
        }, error = function(e) { shiny::showNotification(conditionMessage(e), type = "error", duration = NULL); NULL })
      })
    })
    output$detected <- shiny::renderPlot({ shiny::req(x()); plot_detected_per_sample(x(), input$qc_assay) })
    output$jaccard <- shiny::renderPlot({ shiny::req(x()); x <- x(); ann <- x$annotation[[input$qc_assay]]
      plot_heatmap_matrix(qc_jaccard(x, input$qc_assay), "Jaccard index on detection", groups = stats::setNames(as.character(groups_of(x, input$qc_assay)), ann$sample)) })
    output$pca <- shiny::renderPlot({ shiny::req(x()); plot_pca(x(), input$qc_assay) })
    output$dist <- shiny::renderPlot({ shiny::req(x()); plot_intensity_distributions(x(), input$qc_assay) })
    output$log <- shiny::renderPrint({ shiny::req(x()); cat(x()$log, sep = "\n") })
    x
  })
}
