#' Launch the DiffPTM Shiny application
#'
#' @param ... passed to `shiny::shinyApp` options (e.g. `options = list(port = 3838)`).
#' @export
run_app <- function(...) {
  options(shiny.maxRequestSize = 2 * 1024^3)
  shiny::shinyApp(ui = app_ui(), server = app_server, ...)
}

app_ui <- function() {
  shiny::fluidPage(
    shiny::tags$head(shiny::tags$style(shiny::HTML("
      h1 { color:#000099; } .step-help { color:#555; font-size: 0.9em; }
      .shiny-notification { position: fixed; top: 10px; right: 10px; }"))),
    shiny::titlePanel(shiny::div(shiny::span("DiffPTM", style = "font-weight:600"),
                                 shiny::span(" \u2014 protein-corrected differential PTM analysis", style = "font-weight:300; font-size:0.8em; color:#444")),
                      windowTitle = "DiffPTM"),
    shiny::navlistPanel(id = "step", widths = c(2, 10), well = FALSE,
      shiny::tabPanel("1. Import",     mod_import_ui("import")),
      shiny::tabPanel("2. Design",     mod_design_ui("design")),
      shiny::tabPanel("3. Processing", mod_process_ui("process")),
      shiny::tabPanel("4. Statistics", mod_stats_ui("stats")),
      shiny::tabPanel("5. Results",    mod_results_ui("results")),
      shiny::tabPanel("6. Download",   mod_download_ui("download")),
      shiny::tabPanel("Help",          shiny::includeMarkdown(system.file("app", "help.md", package = "DiffPTM")))
    )
  )
}

app_server <- function(input, output, session) {
  raw     <- mod_import_server("import")                      # list(protein, ptm) assays
  annot   <- mod_design_server("design", raw)                 # reactive ptm_experiment (annotated)
  proc    <- mod_process_server("process", annot)             # reactive processed experiment
  stats   <- mod_stats_server("stats", proc)                  # reactive list(x, res, contrasts)
  mod_results_server("results", stats)
  mod_download_server("download", stats)
}
