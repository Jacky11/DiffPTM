#################################################
# Differential Analysis of Proteomics Data
# Server logic: collects the inputs, calls run_pairwise() (R/pipeline.R) and
# returns a zip file with the results.
#################################################

options(shiny.maxRequestSize = 10 * 1024^3)  # uploads up to 10 GB

# User-friendly error messages
explain_error <- function(msg) {
  if (grepl("read.table", msg)) {
    paste0("Error in read.table : check if the file is comma or tabulation separated. Error message : ", msg)
  } else if (grepl("subscript out of bounds|indice hors limites", msg) || grepl("factor level \\[[0-9]+\\] is duplicated", msg) ||
             grepl("arguments imply differing number of rows", msg)) {
    paste0("Error : check the column names of your protein groups file. Error message : ", msg)
  } else if (grepl("no applicable method for 'ph_with' applied to an object of class ", msg)) {
    paste0("Error : you have to add more Intensities for the two conditions. Error message : ", msg)
  } else msg
}

shinyServer(function(input, output, session) {

  pdf(NULL)  # null graphics device (prevents Rplots.pdf)

  log_messages <- reactiveValues(errors = NULL, stdout = NULL)
  output$error_output <- renderPrint({
    if (!is.null(log_messages$errors)) cat(paste(log_messages$errors, collapse = "\n"))
  })
  output$stdout_output <- renderPrint({
    if (!is.null(log_messages$stdout)) cat(paste(tail(log_messages$stdout, 20), collapse = "\n"))
  })

  #################################
  # Input data
  #################################
  datasetInput <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    tryCatch(
      read_input(inFile$datapath, spectronaut = input$spectronaut, intensity = input$nameKind,
                 name = inFile$name, size = inFile$size),
      error = function(e) { log_messages$errors <- explain_error(conditionMessage(e)); NULL })
  })

  # no iBAQ column in the file: iBAQ values cannot be used
  observeEvent(datasetInput(), {
    if (!any(grepl("^iBAQ\\.", names(datasetInput()$data_init))))
      updateRadioButtons(session, "nameIBAQ", selected = "No")
  })

  #################################
  # Dynamic UI elements
  #################################
  output$protnameDel <- renderUI({
    ds <- datasetInput()
    if (is.null(ds) || ncol(ds$datdelprot) == 0) return(NULL)
    selectInput("varsDelProt", "Columns indicating by + proteins to delete:",
                names(ds$datdelprot), names(ds$datdelprot), multiple = TRUE, width = '800px')
  })

  output$varselectC1 <- renderUI({
    ds <- datasetInput()
    if (is.null(ds)) return(NULL)
    selectInput("vars1", "Choose variables containing intensities for the 1st condition:",
                names(ds$data_fin), names(ds$data_fin), multiple = TRUE)
  })

  output$varselectC2 <- renderUI({
    ds <- datasetInput()
    if (is.null(ds)) return(NULL)
    selectInput("vars2", "Choose variables containing intensities for the 2nd condition:",
                names(ds$data_fin), names(ds$data_fin), multiple = TRUE)
  })

  # Parameters of the analysis, from the inputs
  current_params <- function() {
    utils::modifyList(default_params(), list(
      input_file = input$file1$name,
      spectronaut = isTRUE(input$spectronaut),
      intensity = input$nameKind,
      delete_columns = if (ncol(datasetInput()$datdelprot) == 0) character(0) else as.character(input$varsDelProt),
      cond1_name = input$nameC1, cond1_samples = input$vars1,
      cond2_name = input$nameC2, cond2_samples = input$vars2,
      ibaq = input$nameIBAQ == "Yes",
      min_unique_peptides = input$N_Unique_Identified,
      min_observed = input$numFilter,
      min_percent_observed = input$numFilter2,
      quantified_in = input$quantIn,
      present_min = if (is.null(input$presentMin) || is.na(input$presentMin)) NA else input$presentMin,
      absent_max = if (is.null(input$absentMax) || is.na(input$absentMax)) 0 else input$absentMax,
      normalization = input$nameNorm,
      imputation = input$nameImpute,
      fold_change = input$FCT,
      test = input$nameTest,
      pi0_method = input$nameProp,
      conditional_pvalues = isTRUE(input$condP),
      fdr = input$FDRT,
      saint = isTRUE(input$saint),
      saint_avgp = input$avgp,
      editable_pptx = isTRUE(input$modif),
      overlap = isTRUE(input$overlap),
      overlap_prefix = input$revPrefix,
      overlap_groups = input$groupsManual,
      overlap_nboot = input$nboot,
      overlap_nperm = input$nperm,
      n_cores = input$ncores,
      seed = input$seed
    ), keep.null = TRUE)
  }

  #################################
  # Analysis
  #################################
  result_zip <- reactiveVal(NULL)

  observeEvent(input$run, {
    result_zip(NULL)
    log_messages$errors <- NULL
    log_messages$stdout <- NULL
    ds <- datasetInput()
    if (is.null(ds)) {
      log_messages$errors <- "Please choose an input file first."
      return(NULL)
    }
    last_step <- ""
    out_dir <- file.path(tempdir(), paste0("pairwise_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", sample.int(1e6, 1)))
    res <- withProgress(message = 'Start !', value = 0, {
      tryCatch(
        run_pairwise(current_params(), out_dir = out_dir, data = ds,
                     progress = function(m) incProgress(1 / 16, message = m),
                     log = function(m) {
                       message(m)
                       if (!startsWith(m, "Warning")) last_step <<- m
                       log_messages$stdout <- c(log_messages$stdout, m)
                     }),
        error = function(e) {
          message("Error: ", conditionMessage(e))
          log_messages$errors <- c(log_messages$errors,
                                   paste0(if (nzchar(last_step)) paste0("Error during the step '", last_step, "'. ") else "",
                                          explain_error(conditionMessage(e))))
          NULL
        })
    })
    if (is.null(res)) return(NULL)
    zipfile <- file.path(out_dir, paste0("Differential_Analysis_", input$nameC2, "_vs_", input$nameC1, ".zip"))
    zip::zip(zipfile, files = substring(res$files, nchar(out_dir) + 2), root = out_dir)
    result_zip(zipfile)
  })

  output$download_ui <- renderUI({
    if (is.null(result_zip())) return(NULL)
    tagList(br(), downloadButton("downloadData", "Download Analysis Files"))
  })

  output$downloadData <- downloadHandler(
    filename = function() basename(result_zip()),
    content = function(file) file.copy(result_zip(), file),
    contentType = "application/zip"
  )
})
