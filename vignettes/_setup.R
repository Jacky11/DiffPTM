# shared chunk options of the vignettes
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message = FALSE, warning = FALSE, fig.width = 7, fig.height = 4.5, dpi = 96, out.width = "100%")
fx <- function(...) system.file("extdata", "formats", ..., package = "DiffPTM")
