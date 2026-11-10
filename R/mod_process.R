mod_process_ui <- function(id) {
  ns <- shiny::NS(id)
  help <- function(...) shiny::p(class = "step-help", ...)
  shiny::tagList(
    shiny::h3("Filtering, matching, normalisation"),
    help("The defaults suit most data sets: click ", shiny::strong("Run processing"), " at the bottom right, then check the quality-control plots. ",
         "Every option chosen here is written to the Methods paragraph of the report."),
    shiny::fluidRow(
      shiny::column(4,
        shiny::h4("Filtering"),
        shiny::checkboxInput(ns("flags"), "Remove reverse hits / contaminants", TRUE),
        help("Entries flagged by the software as decoy (reverse), contaminant or 'only identified by site' (MaxQuant '+' columns)."),
        shiny::checkboxInput(ns("merge_redundant"), "Merge site rows repeated for several protein accessions (identical quantities)", TRUE),
        help("Some site reports (e.g. Spectronaut) repeat a site once per protein or isoform with the same values; counting it once avoids inflating the number of tests."),
        shiny::numericInput(ns("min_prob"), "Minimal localization probability", 0.75, 0, 1, 0.05),
        help("Sites whose modified residue is uncertain are removed. 0.75 is the usual 'class I' threshold; ignored for tables without localization probability."),
        shiny::numericInput(ns("min_n"), "Min. observed values in at least one group", 2, 0),
        shiny::numericInput(ns("min_pct"), "Min. % observed values in at least one group", 0, 0, 100),
        help("A feature is kept if it is observed in at least this number (and this percentage) of the samples of at least one group. ",
             "2 is the minimum for a variance; raise it (e.g. 3 of 4) for fewer but more reliable tests.")
      ),
      shiny::column(4,
        shiny::h4("PTM to protein matching"),
        help("For the usage test (DPU), each PTM feature is compared with the abundance of its own protein in the total proteome. ",
             "Matching finds that protein."),
        shiny::radioButtons(ns("match"), NULL, c("By protein identifier (default)" = "id", "By peptide sequence (needs the FASTA)" = "sequence",
                                                  "Hybrid: sequence, then identifier" = "hybrid")),
        shiny::tags$ul(class = "step-help", style = "padding-left: 18px;",
          shiny::tags$li(shiny::strong("Identifier: "), "a site is matched to the protein group(s) of the proteome that share one of its ",
                         "accessions. Enough when both tables come from the same software and database (MaxQuant site tables, ",
                         "Spectronaut and Proteome Discoverer site reports)."),
          shiny::tags$li(shiny::strong("Sequence: "), "the peptide is searched in the protein sequences of the FASTA and the proteins ",
                         "containing it are matched. Use it when the two tables were searched separately and their protein groups ",
                         "or accessions differ (another inference, another database version)."),
          shiny::tags$li(shiny::strong("Hybrid: "), "sequence when it gives a match, identifier otherwise.")),
        shiny::fileInput(ns("fasta"), "FASTA file used for the search (optional)", accept = c(".fasta", ".fa", ".fas", ".txt")),
        help(shiny::strong("When is the FASTA needed? "), "(1) for sequence or hybrid matching; (2) for peptide-level PTM tables that give the ",
             "modified peptide but not the site position (DIA-NN precursor matrix, FragPipe modified peptides): the FASTA locates the ",
             "modified residue in the protein. MaxQuant site tables and site reports of the other software do not need it."),
        shiny::radioButtons(ns("shared"), "Site matching several quantified proteins",
                            c("Average of these proteins" = "mean", "First protein" = "first", "Exclude from the usage test" = "drop")),
        help("A site on a peptide shared by several proteins (or isoforms) can match several protein groups of the proteome. ",
             shiny::strong("Average"), " (default) corrects with the mean of their abundances; ", shiny::strong("First"),
             " uses the first (leading) protein; ", shiny::strong("Exclude"), " keeps such sites for the abundance test (DPA) only."),
        shiny::checkboxInput(ns("sites"), "Collapse peptidoforms to sites (sum of the peptidoforms of a site)", FALSE),
        help("For peptide-level tables: the peptidoforms (charge states, missed cleavages, other modifications) of the same site are summed ",
             "into one feature per site. Leave unticked for site tables (MaxQuant, Spectronaut site reports), which already have one row per site.")
      ),
      shiny::column(4,
        shiny::h4("Normalisation"),
        shiny::radioButtons(ns("norm"), NULL, c("Global median centering (default)" = "global_median", "Per-condition median centering" = "median",
                                                 "ComBat (second factor = batch)" = "combat", "None" = "none")),
        shiny::tags$ul(class = "step-help", style = "padding-left: 18px;",
          shiny::tags$li(shiny::strong("Global median centering: "), "every sample of a table is shifted to the same median. Removes ",
                         "differences of loading and of enrichment efficiency between samples. Use it unless you have a reason not to."),
          shiny::tags$li(shiny::strong("Per-condition median centering: "), "samples are centred within each condition, so a global shift ",
                         "between conditions is kept. Only when a global change of the modification is expected and loading was ",
                         "controlled independently: otherwise a loading difference is reported as a change of every site."),
          shiny::tags$li(shiny::strong("ComBat: "), "removes a batch effect; the batch is the second factor of the Design step."),
          shiny::tags$li(shiny::strong("None: "), "data already normalised by the software (e.g. Spectronaut or DIA-NN normalised quantities; ",
                         "MaxQuant LFQ intensities are normalised between samples but can still be centred).")),
        shiny::h4("Imputation"),
        shiny::radioButtons(ns("impute"), NULL, c("None (recommended)" = "none", "imp4p MLE" = "mle", "imp4p SLSA" = "slsa", "imp4p multiple imputation" = "mi")),
        help("The tests use the observed values. Imputation is optional; it is not needed for the usage test and its effect on the ",
             "results has not been benchmarked here."),
        shiny::actionButton(ns("run"), "Run processing", class = "btn-primary", icon = shiny::icon("gears"))
      )
    ),
    shiny::hr(),
    shiny::h4("Quality control"),
    help("Check that no run stands out (far fewer features, shifted distribution, isolated in the PCA or in the correlation matrix). ",
         "Such runs can be removed at the Design step (empty condition)."),
    shiny::fluidRow(
      shiny::column(3, shiny::selectInput(ns("qc_assay"), "Assay", c("PTM" = "ptm", "Total proteome" = "protein"))),
      shiny::column(9, shiny::verbatimTextOutput(ns("log")))
    ),
    shiny::fluidRow(
      shiny::column(6, shiny::plotOutput(ns("detected"), height = 300), shiny::plotOutput(ns("dist"), height = 300)),
      shiny::column(6, shiny::plotOutput(ns("pca"), height = 380))
    ),
    shiny::fluidRow(
      shiny::column(6, shiny::plotOutput(ns("pearson"), height = 440)),
      shiny::column(6, shiny::plotOutput(ns("jaccard"), height = 440))
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
          if (!is.null(input$fasta)) x$files <- c(x$files, fasta = input$fasta$name)
          x
        }, error = function(e) { shiny::showNotification(conditionMessage(e), type = "error", duration = NULL); NULL })
      })
    })
    output$detected <- shiny::renderPlot({ shiny::req(x()); plot_detected_per_sample(x(), input$qc_assay) })
    output$jaccard <- shiny::renderPlot({ shiny::req(x()); x <- x(); ann <- x$annotation[[input$qc_assay]]
      plot_heatmap_matrix(qc_jaccard(x, input$qc_assay), "Jaccard index on detection", groups = stats::setNames(as.character(groups_of(x, input$qc_assay)), ann$sample)) })
    output$pearson <- shiny::renderPlot({ shiny::req(x()); plot_correlation_matrix(x(), input$qc_assay) })
    output$pca <- shiny::renderPlot({ shiny::req(x()); plot_pca(x(), input$qc_assay) })
    output$dist <- shiny::renderPlot({ shiny::req(x()); plot_intensity_distributions(x(), input$qc_assay) })
    output$log <- shiny::renderPrint({ shiny::req(x()); cat(x()$log, sep = "\n") })
    x
  })
}
