mod_stats_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("Differential analysis"),
    shiny::fluidRow(
      shiny::column(4,
        shiny::checkboxGroupInput(ns("methods"), "Statistical methods (run side by side)",
                                  c("Single model (limma; paired when PTM and proteome samples share names)" = "limma", "MSstatsPTM (model then correct)" = "msstatsptm", "msqrob2PTM (correct then model)" = "msqrob2ptm"),
                                  selected = "limma"),
        shiny::uiOutput(ns("method_note")),
        shiny::p(class = "step-help", "Running time, for 21,770 sites and 4 contrasts on one thread: limma 12 s, msqrob2PTM 1 min, ",
                 "MSstatsPTM 11 min (it fits a model per site and per protein); MSstatsPTM is best kept for a final, confirmatory run."),
        shiny::checkboxGroupInput(ns("analyses"), "Analyses", c("Protein abundance (DPA)" = "protein_DPA", "PTM abundance (DPA)" = "ptm_DPA", "PTM usage, corrected for protein (DPU)" = "ptm_DPU"),
                                  selected = c("protein_DPA", "ptm_DPA", "ptm_DPU")),
        shiny::radioButtons(ns("ctype"), "Contrasts", c("Each condition vs a reference" = "vs_reference", "Each pair of conditions" = "pairwise", "Each condition vs the others" = "vs_all"), inline = TRUE),
        shiny::uiOutput(ns("reference_ui")),
        shiny::radioButtons(ns("blocking"), "Within-sample pairing of site and protein (usage) / subject blocking",
                            c("Per-site pairing, fixed sample effect (default)" = "auto",
                              "Consensus correlation, random sample effect (complete pairs)" = "random", "None (unpaired)" = "none"), inline = FALSE)
      ),
      shiny::column(4,
        shiny::selectInput(ns("pi0"), "pi0 estimator for FDR (cp4p)", c("pounds", "abh", "st.spline", "st.boot", "langaas", "histo", "jiang", "slim", "bh")),
        shiny::numericInput(ns("fdr"), "FDR threshold", 0.05, 0, 1, 0.01),
        shiny::numericInput(ns("fc"), "Minimal fold change (linear; 1 = no fold-change filter)", 1, 1, step = 0.1),
        shiny::actionButton(ns("run"), "Run", class = "btn-primary", icon = shiny::icon("play")),
        shiny::uiOutput(ns("status"))
      ),
      shiny::column(4, shiny::h4("Contrasts to be tested"),
                    shiny::p(class = "step-help", "With a second factor (e.g. time), conditions are compared within each of its levels."),
                    shiny::tableOutput(ns("contrasts")))
    )
  )
}

mod_stats_server <- function(id, processed) {
  shiny::moduleServer(id, function(input, output, session) {
    output$method_note <- shiny::renderUI({
      miss <- setdiff(c("msstatsptm", "msqrob2ptm"), available_methods())
      if (length(miss)) shiny::p(class = "step-help", "Not installed on this server: ", paste(miss, collapse = ", "))
    })
    output$reference_ui <- shiny::renderUI({
      x <- processed(); shiny::req(x, input$ctype == "vs_reference")
      cc <- unique(x$annotation$ptm$condition)
      ctl <- grep("^(mock|ctrl|control|ctl|cont|wt|wildtype|untreated|untr|vehicle|veh|dmso|ns|unstim|neg|basal|baseline|ref|t0|0h)$", cc, ignore.case = TRUE, value = TRUE)
      shiny::selectInput(session$ns("reference"), "Reference (control) condition", cc, selected = if (length(ctl)) ctl[1] else cc[1])
    })
    contrasts <- shiny::reactive({ shiny::req(processed())
      tryCatch(make_contrasts(processed(), input$ctype, reference = if (input$ctype == "vs_reference") input$reference else NULL),
               error = function(e) NULL) })
    output$contrasts <- shiny::renderTable(contrasts())
    out <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$run, {               # runs immediately (an eventReactive would wait for the Results page)
      x <- processed(); shiny::req(x)
      sel <- intersect(input$methods, available_methods())
      shiny::withProgress(message = "Fitting models", {
        shiny::req(contrasts())
        res <- tryCatch(run_methods(x, contrasts(), methods = sel, analyses = input$analyses, pi0_method = input$pi0, blocking = input$blocking),
                        error = function(e) { shiny::showNotification(conditionMessage(e), type = "error", duration = NULL); NULL })
      })
      shiny::req(res)
      out(list(x = x, res = res, contrasts = contrasts(), fdr = input$fdr, min_log2fc = log2(input$fc), methods = sel, blocking = input$blocking))
    })
    output$status <- shiny::renderUI({
      s <- out(); shiny::req(s)
      r <- s$res; r$sig <- !is.na(r$adj_pvalue) & r$adj_pvalue < s$fdr & abs(r$log2FC) >= s$min_log2fc
      tab <- stats::aggregate(sig ~ analysis + method, r, sum)
      lab <- c(protein_DPA = "proteins (DPA)", ptm_DPA = "PTM features (DPA)", ptm_DPU = "PTM features, usage (DPU)")
      shiny::div(class = "alert alert-success", style = "padding: 8px 12px; margin-top: 10px;",
        shiny::icon("circle-check"), shiny::strong(" Analysis done. "),
        sprintf("Significant at FDR %s, all contrasts: ", s$fdr),
        shiny::tags$ul(style = "margin-bottom: 4px;", lapply(seq_len(nrow(tab)), function(i)
          shiny::tags$li(sprintf("%s, %s: %d", tab$method[i], lab[[tab$analysis[i]]], tab$sig[i])))),
        shiny::strong("See 5. Results."))
    })
    out
  })
}
