#################################################
# Installs the R packages needed by the application
#   source("scripts/install_packages.R")
#################################################

cran <- c("shiny", "officer", "openxlsx", "yaml", "zip",
          "imp4p", "cp4p", "statmod",
          "ggplot2", "ggrepel", "ggridges", "ggdendro", "ggplotify", "gridExtra", "reshape2", "dplyr",
          "VennDiagram", "venn", "UpSetR", "pvclust", "vegan",
          "export")                               # optional: editable PowerPoint plots
bioc <- c("limma", "multtest", "qvalue")          # limma; multtest and qvalue are needed by cp4p

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, update = FALSE)
missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) install.packages(missing_cran)
