#################################################
# Loads the packages and the functions of the application.
# Used by the Shiny app (global.R) and by scripts/run_analysis.R:
#   source("load_app.R")
#################################################

local({
  required <- c(
    "shiny", "officer", "openxlsx", "yaml", "zip",          # app, outputs
    "imp4p", "cp4p", "limma", "statmod",                    # statistics
    "ggplot2", "ggrepel", "ggridges", "ggdendro", "ggplotify", "gridExtra", "reshape2", "dplyr",
    "VennDiagram", "venn", "UpSetR",                        # R/venn.R
    "pvclust", "vegan"                                      # R/overlap_clustering.R
  )
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0)
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         "\nInstall them with: source('scripts/install_packages.R')", call. = FALSE)
})

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(ggridges)
  library(officer)
  library(limma)
  library(statmod)
})

local({
  # folder of the app = folder of this file
  app_dir <- getOption("pairwise.app_dir")
  if (is.null(app_dir)) {
    f <- NULL
    for (i in rev(seq_len(sys.nframe()))) {
      f <- sys.frame(i)$ofile
      if (!is.null(f)) break
    }
    app_dir <- if (is.null(f)) getwd() else dirname(normalizePath(f))
    options(pairwise.app_dir = app_dir)
  }
  for (f in c("adjustp.R", "limma-test.R", "venn.R", "overlap_clustering.R", "utils.R", "params.R", "filters.R", "read_data.R", "saint.R", "pipeline.R"))
    sys.source(file.path(app_dir, "R", f), envir = globalenv())
})
