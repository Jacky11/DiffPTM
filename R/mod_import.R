mod_import_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    # the browser knows the modification date of the chosen files; send it with the name and size (for the Methods paragraph)
    shiny::tags$script(shiny::HTML(sprintf(
      "$(document).on('change', '#%s, #%s', function(e) { var f = e.target.files && e.target.files[0];
         if (f) Shiny.setInputValue(e.target.id + '_meta', {name: f.name, size: f.size, modified: f.lastModified}); });",
      ns("protein_file"), ns("ptm_file")))),
    shiny::h3("Input tables"),
    shiny::p(class = "step-help", "Upload the total-proteome table and the PTM table produced by your quantification software (MaxQuant, Spectronaut, Proteome Discoverer, FragPipe, DIA-NN), or load the example dataset."),
    shiny::fluidRow(
      shiny::column(4,
        shiny::selectInput(ns("software"), "Software", names(column_specs())),
        shiny::uiOutput(ns("software_help")),
        shiny::fileInput(ns("protein_file"), "Total proteome table (proteins)", accept = c(".txt", ".tsv", ".csv", ".xlsx")),
        shiny::fileInput(ns("ptm_file"), "PTM table (sites or modified peptides)", accept = c(".txt", ".tsv", ".csv", ".xlsx")),
        shiny::checkboxInput(ns("log2"), "Intensities are raw (apply log2)", TRUE),
        shiny::conditionalPanel(sprintf("input['%s'] == 'MaxQuant'", ns("software")),
          shiny::radioButtons(ns("multiplicity"), "Multiplicity columns (___1/2/3)",
                              c("One feature per multiplicity" = "expand", "Keep ___1 only" = "first", "Sum" = "sum"), inline = TRUE)),
        shiny::actionButton(ns("example"), "Load example dataset", icon = shiny::icon("flask")),
        shiny::actionButton(ns("go"), "Import", class = "btn-primary", icon = shiny::icon("upload"))
      ),
      shiny::column(8,
        shiny::uiOutput(ns("status")),
        shiny::verbatimTextOutput(ns("summary")),
        shiny::tags$details(
          shiny::tags$summary(shiny::strong("Advanced: column mapping"), " (only if the quantitative columns are not found)"),
          shiny::p(class = "step-help", "Regular expressions with one capture group giving the sample name. Leave empty to use the defaults of the selected software."),
          shiny::textInput(ns("protein_quant"), "Protein quantitative columns (regex)", ""),
          shiny::textInput(ns("ptm_quant"), "PTM quantitative columns (regex)", ""))
      )
    )
  )
}

mod_import_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    assays <- shiny::reactiveVal(NULL)
    failed <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$go, {
      shiny::req(input$protein_file, input$ptm_file)
      shiny::withProgress(message = "Importing", {
        res <- tryCatch({
          qp <- if (nzchar(input$protein_quant)) input$protein_quant else NULL
          qq <- if (nzchar(input$ptm_quant)) input$ptm_quant else NULL
          list(
            protein = import_protein(input$protein_file$datapath, input$software, quant_pattern = qp, log2_transform = input$log2),
            ptm = import_ptm(input$ptm_file$datapath, input$software, quant_pattern = qq, log2_transform = input$log2,
                             multiplicity = input$multiplicity %||% "expand"))
        }, error = function(e) { failed(conditionMessage(e)); shiny::showNotification(conditionMessage(e), type = "error", duration = NULL); NULL })
      })
      if (!is.null(res)) {
        failed(NULL)
        for (k in c("protein", "ptm")) {                  # name, size and date of the file on the user's computer
          fi <- input[[paste0(k, "_file")]]; mt <- input[[paste0(k, "_file_meta")]]
          mod <- if (!is.null(mt) && identical(mt$name, fi$name) && !is.null(mt$modified)) as.POSIXct(mt$modified / 1000, origin = "1970-01-01") else NA
          res[[k]]$source <- file_source(fi$datapath, input$software, name = fi$name, size = fi$size, modified = mod)
        }
      }
      assays(res)
    })
    shiny::observeEvent(input$example, {
      ex <- system.file("extdata", "example", package = "DiffPTM")
      if (!nzchar(ex) || !file.exists(file.path(ex, "protein.tsv"))) {
        shiny::showNotification("Example data not installed; see inst/extdata/example/README.md", type = "warning"); return()
      }
      assays(list(protein = import_protein(file.path(ex, "protein.tsv"), "Generic"),
                  ptm = import_ptm(file.path(ex, "ptm.tsv"), "Generic"),
                  annotation = read_annotation(file.path(ex, "annotation.tsv"))))
      failed(NULL)
      shiny::showNotification("Example dataset loaded (simulated, 2 conditions x 4 samples); go to step 2.", type = "message")
    })
    output$software_help <- shiny::renderUI({
      txt <- switch(input$software,
        MaxQuant = "proteinGroups.txt and Phospho (STY)Sites.txt / GlyGly (K)Sites.txt. LFQ intensities are used when present. In a combined search all runs appear in every table; step 2 keeps the right ones and pairs them automatically.",
        Spectronaut = "Protein pivot report (PG.Quantity) and PTM site pivot report (PTM.Quantity, PTM.SiteProbability); semicolon-separated files with decimal commas are read as such.",
        `Proteome Discoverer` = "Proteins export and PeptideGroups export (tab-separated text); sites are read from 'Modifications in Master Proteins'.",
        FragPipe = "combined_protein.tsv or abundance_protein_*.tsv, and abundance_single-site_*.tsv (TMT) or combined_modified_peptide.tsv (label-free).",
        `DIA-NN` = "report.pg_matrix.tsv and report.pr_matrix.tsv (or a site matrix); sample names are the raw file names.",
        Generic = "Wide tables: protein_id / feature_id columns, metadata columns and one intensity_<sample> column per sample (see Help).")
      shiny::p(class = "step-help", txt)
    })
    output$status <- shiny::renderUI({
      a <- assays(); err <- failed()
      if (is.null(a) && is.null(err)) return(shiny::div(class = "alert alert-info", style = "padding: 8px 12px;",
        "No data imported yet: choose the two files and click Import (or load the example dataset)."))
      if (is.null(a)) return(shiny::div(class = "alert alert-danger", style = "padding: 8px 12px;",
        shiny::icon("triangle-exclamation"), shiny::strong(" Import failed. "), err))
      n <- function(k) sprintf("%d %s x %d runs", nrow(a[[k]]$intensity), if (k == "protein") "proteins" else "PTM features", ncol(a[[k]]$intensity))
      src <- a$ptm$source
      shiny::div(class = "alert alert-success", style = "padding: 8px 12px;",
        shiny::icon("circle-check"), shiny::strong(" Both tables imported successfully"),
        sprintf("%s. Proteome: %s; PTM: %s. ", if (!is.null(src)) sprintf(" (%s: %s and %s)", src$software, a$protein$source$name, src$name) else "",
                n("protein"), n("ptm")),
        shiny::strong("Next: 2. Design."))
    })
    output$summary <- shiny::renderPrint({
      a <- assays(); shiny::req(a)
      for (n in intersect(names(a), c("protein", "ptm"))) cat(sprintf("%s: %d features x %d samples\n   samples: %s\n", n,
                                      nrow(a[[n]]$intensity), ncol(a[[n]]$intensity),
                                      paste(utils::head(colnames(a[[n]]$intensity), 8), collapse = ", ")))
    })
    assays
  })
}
