mod_design_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("Experimental design"),
    shiny::p(class = "step-help",
             "Assign each sample to a condition, an optional second factor (time, batch) and an optional subject for paired designs. ",
             "Samples left without a condition are excluded. When the PTM and proteome samples of the same biological sample carry ",
             "different names (e.g. Phospho_mock_1 and mock_1), give them the same name with the patterns below: the DPU test then pairs them."),
    shiny::fluidRow(
      shiny::column(4,
        shiny::h4("1. Sample names"),
        shiny::p(class = "step-help", "Optional: a regular expression per table selects the samples (those that match) and names them (the capture group). Leave empty to keep all samples as they are."),
        shiny::textInput(ns("rename_protein"), "Proteome samples: regex, capture group = sample name", "", placeholder = "e.g. ^Proteome-\\d+_(.+)$"),
        shiny::textInput(ns("rename_ptm"), "PTM samples: regex, capture group = sample name", "", placeholder = "e.g. ^Phospho_(.+)$"),
        shiny::verbatimTextOutput(ns("rename_preview")),
        shiny::h4("2. Conditions"),
        shiny::textInput(ns("pattern"), "Token separator in sample names (regex)", "[_ .-]"),
        shiny::textInput(ns("cond_pos"), "Token(s) used as condition (-1 = all but the last; or 1, 1-2, ...)", "-1"),
        shiny::textInput(ns("factor_pos"), "Token(s) used as second factor (empty = none)", ""),
        shiny::textInput(ns("subject_pos"), "Token(s) used as subject (empty = none)", ""),
        shiny::actionButton(ns("guess"), "Pre-fill from sample names"),
        shiny::fileInput(ns("ann_file"), "...or upload an annotation table (tsv/csv: sample, condition, factor, subject, mixture)"),
        shiny::checkboxInput(ns("same"), "Proteome and PTM samples are the same (copy the proteome table)", FALSE),
        shiny::actionButton(ns("apply"), "Apply design", class = "btn-primary")
      ),
      shiny::column(8,
        shiny::p(class = "step-help", "Edit cells by double-clicking. Empty condition = sample excluded."),
        shiny::h4("Total proteome samples"), DT::DTOutput(ns("ann_protein")),
        shiny::h4("PTM samples"), DT::DTOutput(ns("ann_ptm")),
        shiny::verbatimTextOutput(ns("status"))
      )
    )
  )
}

mod_design_server <- function(id, raw) {
  shiny::moduleServer(id, function(input, output, session) {
    ann <- shiny::reactiveValues(protein = NULL, ptm = NULL)
    blank <- function(samples) data.frame(sample = samples, condition = "", factor = "", subject = "", mixture = "", stringsAsFactors = FALSE)
    pat <- function(k) { p <- input[[paste0("rename_", k)]]; if (is.null(p) || !nzchar(p)) "" else p }
    # assays with the renaming applied (the tables below refer to the new names)
    renamed <- shiny::reactive({
      a <- raw(); shiny::req(a)
      for (k in c("protein", "ptm")) a[[k]] <- tryCatch(rename_samples(a[[k]], pat(k)), error = function(e) a[[k]])
      a
    })
    output$rename_preview <- shiny::renderText({
      a <- raw(); shiny::req(a); r <- renamed()
      lines <- vapply(c("protein", "ptm"), function(k) {
        old <- colnames(a[[k]]$intensity); new <- colnames(r[[k]]$intensity)
        chk <- tryCatch({ rename_samples(a[[k]], pat(k)); "" }, error = function(e) paste(" !!", conditionMessage(e)))
        sprintf("%s: %d of %d samples kept, e.g. %s -> %s%s", k, length(new), length(old), old[1], new[1], chk) }, "")
      paste(c(lines, sprintf("%d sample names shared by both tables (paired analysis possible)", length(intersect(colnames(r$protein$intensity), colnames(r$ptm$intensity))))), collapse = "\n")
    })
    reset_tables <- function() { r <- renamed(); ann$protein <- blank(colnames(r$protein$intensity)); ann$ptm <- blank(colnames(r$ptm$intensity)) }
    shiny::observeEvent(raw(), {
      a <- raw(); shiny::req(a)
      reset_tables()
      if (!is.null(a$annotation)) for (k in c("protein", "ptm")) {       # example dataset comes with its annotation
        d <- a$annotation; s <- colnames(a[[k]]$intensity); m <- d[match(s, d$sample), ]
        m[is.na(m)] <- ""; m$sample <- s; ann[[k]] <- m[, names(blank(s))]
      }
    })
    shiny::observeEvent(c(input$rename_protein, input$rename_ptm), { shiny::req(raw()); reset_tables() }, ignoreInit = TRUE)
    shiny::observeEvent(input$guess, {
      r <- renamed(); shiny::req(r)
      for (k in c("protein", "ptm")) {
        g <- guess_annotation(colnames(r[[k]]$intensity), input$pattern, input$cond_pos, input$factor_pos, input$subject_pos)
        g[is.na(g)] <- ""; ann[[k]] <- g
      }
    })
    shiny::observeEvent(input$ann_file, {
      d <- read_annotation(input$ann_file$datapath); r <- renamed()
      for (k in c("protein", "ptm")) { s <- colnames(r[[k]]$intensity); m <- d[match(s, d$sample), ]; m[is.na(m)] <- ""; m$sample <- s; ann[[k]] <- m }
    })
    output$ann_protein <- DT::renderDT(ann$protein, editable = TRUE, options = list(pageLength = 12, dom = "tp"), rownames = FALSE)
    output$ann_ptm     <- DT::renderDT(ann$ptm, editable = TRUE, options = list(pageLength = 12, dom = "tp"), rownames = FALSE)
    shiny::observeEvent(input$ann_protein_cell_edit, ann$protein <- DT::editData(ann$protein, input$ann_protein_cell_edit, rownames = FALSE))
    shiny::observeEvent(input$ann_ptm_cell_edit,     ann$ptm     <- DT::editData(ann$ptm,     input$ann_ptm_cell_edit,     rownames = FALSE))

    x <- shiny::eventReactive(input$apply, {
      r <- renamed(); shiny::req(r)
      if (input$same) { ann$ptm <- ann$protein[ann$protein$sample %in% colnames(r$ptm$intensity), ] }
      tryCatch({
        mk <- function(d) {
          d <- d[!is.na(d$condition) & nzchar(trimws(d$condition)), , drop = FALSE]
          if (!nrow(d)) stop("No sample has a condition")
          d[d == ""] <- NA
          make_annotation(d$sample, d$condition, d$factor, d$subject, d$mixture)
        }
        ap <- mk(ann$protein); aq <- mk(ann$ptm)
        new_ptm_experiment(keep_samples(r$protein, ap$sample), keep_samples(r$ptm, aq$sample),
                           annotation = list(protein = ap, ptm = aq))
      }, error = function(e) { shiny::showNotification(conditionMessage(e), type = "error", duration = NULL); NULL })
    })
    output$status <- shiny::renderPrint({ shiny::req(x()); print(x()) })
    x
  })
}
