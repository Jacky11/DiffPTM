mod_design_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("Experimental design"),
    shiny::p(class = "step-help",
             "DiffPTM proposes a design from the imported tables. Check it, correct it in the tables if needed, then click ",
             shiny::strong("Apply design"), "."),
    shiny::fluidRow(
      shiny::column(5,
        shiny::h4("1. What DiffPTM found"),
        shiny::uiOutput(ns("suggestion")),
        shiny::h4("2. Conditions"),
        shiny::p(class = "step-help",
                 "Number of samples per condition in each table, and how many PTM samples are paired with a proteome ",
                 "sample (same ", shiny::em("sample"), " name)."),
        shiny::tableOutput(ns("summary")),
        shiny::h4("3. Apply"),
        shiny::actionButton(ns("apply"), "Apply design", class = "btn-primary"),
        shiny::verbatimTextOutput(ns("status")),
        shiny::tags$details(
          shiny::tags$summary(shiny::strong("Advanced options"), " (rename runs in bulk, read conditions differently, annotation file)"),
          shiny::div(style = "margin-top: 10px;",
            shiny::p(class = "step-help",
                     "Rename in bulk: a regular expression applied to the run names; its capture group becomes the sample ",
                     "name (e.g. ^P_(.+)$ turns P_R1-A into R1-A). Runs that do not match keep their current sample name."),
            shiny::fluidRow(
              shiny::column(6, shiny::textInput(ns("rename_protein"), "Proteome runs", "", placeholder = "e.g. ^T_(.+)$")),
              shiny::column(6, shiny::textInput(ns("rename_ptm"), "PTM runs", "", placeholder = "e.g. ^P_(.+)$"))),
            shiny::actionButton(ns("do_rename"), "Rename"),
            shiny::hr(),
            shiny::p(class = "step-help",
                     "Conditions from sample names: names are split into parts at the separator; choose the parts that form ",
                     "the condition (-1 = all but the last part; 1 = first part; 1-2 = first two parts; 2,3 = parts 2 and 3), ",
                     "and optionally a second factor (e.g. time) and a subject."),
            shiny::fluidRow(
              shiny::column(6, shiny::textInput(ns("pattern"), "Separator (regular expression)", "[_ .-]")),
              shiny::column(6, shiny::textInput(ns("cond_pos"), "Parts used as condition", "-1"))),
            shiny::fluidRow(
              shiny::column(6, shiny::textInput(ns("factor_pos"), "Second factor (empty = none)", "")),
              shiny::column(6, shiny::textInput(ns("subject_pos"), "Subject (empty = none)", ""))),
            shiny::actionButton(ns("guess"), "Read conditions from the sample names"),
            shiny::hr(),
            shiny::fileInput(ns("ann_file"), "Upload an annotation table (tsv/csv; columns run or sample, condition, optional factor, subject, mixture)"),
            shiny::actionButton(ns("auto"), "Restore the automatic suggestion")
          )
        )
      ),
      shiny::column(7,
        shiny::p(class = "step-help",
                 shiny::strong("run"), " = column of the imported table; ", shiny::strong("sample"),
                 " = name used in the analysis (a PTM and a proteome run with the same sample name are the same biological ",
                 "sample and are paired); ", shiny::strong("condition"), " empty = run not used. Double-click a cell to edit it."),
        shiny::h4("Proteome runs"), DT::DTOutput(ns("ann_protein")),
        shiny::h4("PTM runs"), DT::DTOutput(ns("ann_ptm")),
        shiny::p(class = "step-help", style = "margin-top: 8px;",
                 shiny::strong("Optional columns: "), "factor = a second factor such as time or batch (contrasts are then built ",
                 "within each level); subject = the same individual or culture across conditions (paired design); ",
                 "mixture = TMT plex.")
      )
    )
  )
}

mod_design_server <- function(id, raw) {
  shiny::moduleServer(id, function(input, output, session) {
    ann <- shiny::reactiveValues(protein = NULL, ptm = NULL, sugg = NULL)
    cols <- c("run", "sample", "condition", "factor", "subject", "mixture")
    clean <- function(d) { for (cc in setdiff(cols, names(d))) d[[cc]] <- ""; d <- d[, cols]; d[is.na(d)] <- ""; rownames(d) <- NULL; d }
    used <- function(d) d[nzchar(trimws(d$condition)) & nzchar(trimws(d$sample)), , drop = FALSE]

    restore <- function() { s <- ann$sugg; for (k in c("protein", "ptm")) ann[[k]] <- clean(s$annotation[[k]]) }
    shiny::observeEvent(raw(), {
      a <- raw(); shiny::req(a)
      ann$sugg <- suggest_design(a$protein, a$ptm)
      restore()
      if (!is.null(a$annotation)) for (k in c("protein", "ptm")) {       # the example dataset comes with its annotation
        d <- ann[[k]]; m <- match(d$run, a$annotation$sample)
        for (cc in c("condition", "factor", "subject", "mixture")) if (cc %in% names(a$annotation)) {
          v <- as.character(a$annotation[[cc]][m]); d[[cc]][!is.na(v)] <- v[!is.na(v)] }
        ann[[k]] <- clean(d)
      }
    })
    shiny::observeEvent(input$auto, { shiny::req(ann$sugg); restore() })
    shiny::observeEvent(input$do_rename, {
      shiny::req(ann$protein)
      for (k in c("protein", "ptm")) {
        p <- input[[paste0("rename_", k)]]
        if (is.null(p) || !nzchar(p)) next
        d <- ann[[k]]
        hit <- tryCatch(grepl(p, d$run, perl = TRUE), error = function(e) { shiny::showNotification(conditionMessage(e), type = "error"); NULL })
        if (is.null(hit)) next
        if (!any(hit)) { shiny::showNotification(sprintf("No %s run matches %s", k, p), type = "warning"); next }
        d$sample[hit] <- sub(p, "\\1", d$run[hit], perl = TRUE)
        ann[[k]] <- d
      }
    })
    shiny::observeEvent(input$guess, {
      shiny::req(ann$protein)
      for (k in c("protein", "ptm")) {
        d <- ann[[k]]; u <- nzchar(d$sample)
        if (!any(u)) next
        g <- tryCatch(suppressWarnings(guess_annotation(d$sample[u], input$pattern, input$cond_pos, input$factor_pos, input$subject_pos)),
                      error = function(e) NULL)
        if (is.null(g)) next
        for (cc in c("condition", "factor", "subject")) { v <- g[[cc]]; v[is.na(v)] <- ""; d[[cc]][u] <- v }
        ann[[k]] <- d
      }
    })
    shiny::observeEvent(input$ann_file, {
      f <- data.table::fread(input$ann_file$datapath, data.table = FALSE, colClasses = "character")
      names(f) <- tolower(names(f))
      if (!"condition" %in% names(f) || !any(c("run", "sample") %in% names(f))) {
        shiny::showNotification("The annotation table needs a 'run' or 'sample' column and a 'condition' column.", type = "error"); return()
      }
      n_hit <- 0
      for (k in c("protein", "ptm")) {
        d <- ann[[k]]
        m <- if ("run" %in% names(f)) match(d$run, f$run) else match(d$sample, f$sample)
        if ("run" %in% names(f) && all(is.na(m)) && "sample" %in% names(f)) m <- match(d$sample, f$sample)
        n_hit <- n_hit + sum(!is.na(m))
        if ("run" %in% names(f) && "sample" %in% names(f)) d$sample[!is.na(m)] <- f$sample[m[!is.na(m)]]
        for (cc in c("condition", "factor", "subject", "mixture")) {
          v <- if (cc %in% names(f)) f[[cc]][m] else rep(NA, nrow(d)); v[is.na(v)] <- ""
          d[[cc]] <- ifelse(is.na(m), "", v)
        }
        ann[[k]] <- d
      }
      shiny::showNotification(sprintf("%d runs found in the annotation table; the others are not used.", n_hit), type = "message")
    })

    output$suggestion <- shiny::renderUI({
      s <- ann$sugg; shiny::req(s)
      shiny::wellPanel(style = "padding: 10px 15px;",
        shiny::tags$ul(style = "padding-left: 18px; margin-bottom: 6px;", lapply(s$messages, shiny::tags$li)),
        shiny::p(class = "step-help", style = "margin-bottom: 0;",
                 "Paired samples are taken to be the enriched and total fractions of the same biological sample. If they are not, ",
                 "choose 'None (unpaired)' at the Statistics step."))
    })
    output$summary <- shiny::renderTable({
      shiny::req(ann$protein, ann$ptm)
      up <- used(ann$protein); uq <- used(ann$ptm)
      cc <- sort(unique(c(up$condition, uq$condition)))
      if (!length(cc)) return(data.frame(Message = "No condition yet: fill the condition column, or use Advanced options."))
      paired <- intersect(up$sample, uq$sample)
      data.frame(Condition = cc,
                 `Proteome samples` = vapply(cc, function(x) sum(up$condition == x), 0L),
                 `PTM samples` = vapply(cc, function(x) sum(uq$condition == x), 0L),
                 `Paired` = vapply(cc, function(x) sum(uq$sample %in% paired & uq$condition == x), 0L),
                 check.names = FALSE)
    }, digits = 0)
    dt <- function(d) DT::datatable(d, editable = list(target = "cell", disable = list(columns = 0)), rownames = FALSE,
                                    options = list(dom = "t", paging = FALSE, scrollY = "330px", scrollCollapse = TRUE, ordering = FALSE))
    output$ann_protein <- DT::renderDT(dt(ann$protein))
    output$ann_ptm     <- DT::renderDT(dt(ann$ptm))
    shiny::observeEvent(input$ann_protein_cell_edit, ann$protein <- DT::editData(ann$protein, input$ann_protein_cell_edit, rownames = FALSE))
    shiny::observeEvent(input$ann_ptm_cell_edit,     ann$ptm     <- DT::editData(ann$ptm,     input$ann_ptm_cell_edit,     rownames = FALSE))

    x <- shiny::eventReactive(input$apply, {
      a <- raw(); shiny::req(a, ann$protein)
      tryCatch({
        mk <- function(d, assay, what) {
          d <- used(d)
          if (!nrow(d)) stop("No ", what, " run has a sample name and a condition")
          if (anyDuplicated(d$sample)) stop("Duplicated ", what, " sample names: ", paste(unique(d$sample[duplicated(d$sample)]), collapse = ", "))
          d[d == ""] <- NA
          list(assay = set_run_names(keep_samples(assay, d$run), stats::setNames(d$sample, d$run)),
               ann = make_annotation(d$sample, d$condition, d$factor, d$subject, d$mixture))
        }
        p <- mk(ann$protein, a$protein, "proteome"); q <- mk(ann$ptm, a$ptm, "PTM")
        new_ptm_experiment(p$assay, q$assay, annotation = list(protein = p$ann, ptm = q$ann))
      }, error = function(e) { shiny::showNotification(conditionMessage(e), type = "error", duration = NULL); NULL })
    })
    output$status <- shiny::renderPrint({
      e <- x(); shiny::req(e); an <- e$annotation
      cat(sprintf("Design applied: %d proteome and %d PTM samples, %d conditions, %d paired samples.\nGo to 3. Processing.",
                  nrow(an$protein), nrow(an$ptm), length(unique(c(an$protein$condition, an$ptm$condition))),
                  length(intersect(an$protein$sample, an$ptm$sample))))
    })
    x
  })
}
