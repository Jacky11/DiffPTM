#!/usr/bin/env Rscript
#################################################
# Runs an analysis without the Shiny interface, from a YAML parameter file.
#
#   Rscript scripts/run_analysis.R params/example.yml [output_folder]
#
# The parameter file can be written by hand (see params/example.yml) or taken
# from the parameters.yml file saved with the results of any analysis made with
# the Shiny application.
#
# Helpers to fill a parameter file:
#   Rscript scripts/run_analysis.R --list-samples proteinGroups.txt [Raw|LFQ]
#################################################

args <- commandArgs(trailingOnly = TRUE)
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
app_dir <- if (length(this_file) == 1) normalizePath(file.path(dirname(this_file), "..")) else getwd()
options(pairwise.app_dir = app_dir)
source(file.path(app_dir, "load_app.R"))

if (length(args) == 0) {
  cat("Usage: Rscript scripts/run_analysis.R <parameters.yml> [output_folder]\n",
      "       Rscript scripts/run_analysis.R --list-samples <file> [Raw|LFQ] [--spectronaut]\n")
  quit(status = 1)
}

if (args[1] == "--list-samples") {
  intensity <- if (length(args) >= 3 && args[3] %in% c("Raw", "LFQ")) args[3] else "Raw"
  ds <- read_input(args[2], spectronaut = "--spectronaut" %in% args, intensity = intensity)
  cat("Sample columns (", intensity, "):\n  ", paste(names(ds$data_fin), collapse = "\n   "), "\n", sep = "")
  cat("Columns flagging proteins to delete with '+':\n  ", paste(names(ds$datdelprot), collapse = "\n   "), "\n", sep = "")
  quit(status = 0)
}

p <- read_params(args[1])
out_dir <- if (length(args) >= 2) args[2] else
  file.path("results", paste0(p$cond2_name, "_vs_", p$cond1_name))

res <- run_pairwise(p, out_dir = out_dir)
cat("\nNumber of proteins:\n"); print(res$counts)
cat("\nFiles written in", normalizePath(out_dir), ":\n ", paste(basename(res$files), collapse = "\n  "), "\n")
