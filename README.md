# DiffPTM

<!-- badges: start -->
[![R-CMD-check](https://github.com/Jacky11/DiffPTM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Jacky11/DiffPTM/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](LICENSE.md)
<!-- badges: end -->

**Differential analysis of post-translational modifications (PTMs), corrected for protein abundance.**

A modified site can change because the modification changes or because its protein changes.
DiffPTM separates the two:

* **differential PTM abundance (DPA)**: does the site intensity change?
* **differential PTM usage (DPU)**: does the site change *relative to its protein*?

It runs three statistical model forms for the usage test on the same processed data —
a single linear model of site and protein (limma), model-then-correct (MSstatsPTM) and
correct-then-model (msqrob2PTM) — and reports diagnostics that tell which one your data
support. Everything is available from a graphical interface that runs on your own computer
(no data leave it) and from R functions.

> **Development version (0.0.1).** DiffPTM is under active development and the accompanying
> article is in preparation. Functions and defaults may still change; see [NEWS.md](NEWS.md).

## Installation

DiffPTM needs R ≥ 4.2. Some dependencies come from Bioconductor, so install through
`BiocManager`:

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install("Jacky11/DiffPTM", update = FALSE)
```

Install in a fresh R session (in RStudio: *Session → Restart R*): if DiffPTM is already
loaded, R reports that the package "is in use" and does not replace it. `update = FALSE` skips
the offer to update every other installed package; check afterwards with
`packageVersion("DiffPTM")`.

This installs DiffPTM and its required dependencies. The optional packages that add the
MSstatsPTM and msqrob2PTM model forms, imputation, sequence-based site-to-protein matching and
the PowerPoint export can be installed with:

```r
BiocManager::install(c("MSstatsPTM", "msqrob2", "QFeatures", "SummarizedExperiment",
                       "sva", "Biostrings"))
install.packages(c("imp4p", "cp4p", "AhoCorasickTrie", "rvg", "plotly", "ggrepel",
                   "ggdendro", "ggridges", "readxl", "markdown"))
```

To install with the vignettes (slower):

```r
BiocManager::install("Jacky11/DiffPTM", build_vignettes = TRUE, dependencies = TRUE, update = FALSE)
```

Without R: a Docker image serves the interface on <http://localhost:3838>
(build it from this repository: `docker build -t diffptm . && docker run -p 3838:3838 diffptm`).

## Quick start: the interface

```r
DiffPTM::run_app()
```

The interface opens in your browser and guides you through six steps: import (MaxQuant,
Spectronaut, Proteome Discoverer, FragPipe, DIA-NN or a generic table), sample annotation,
processing (filtering, site-to-protein matching, normalisation), quality control, statistics
(DPA and DPU, with any of the three model forms) and export (Excel workbook with an
automatically written Methods paragraph, PowerPoint and PDF). The **Load example dataset**
button on the first page loads a small simulated experiment to try every step.

## Quick start: R

```r
library(DiffPTM)
ex <- system.file("extdata", "example", package = "DiffPTM")
prot <- import_protein(file.path(ex, "protein.tsv"), "Generic")
ptm  <- import_ptm(file.path(ex, "ptm.tsv"), "Generic")
ann  <- read_annotation(file.path(ex, "annotation.tsv"))   # columns: sample, condition

x <- new_ptm_experiment(prot, ptm, annotation = list(protein = ann, ptm = ann)) |>
  filter_localization(0.75) |>
  filter_min_observed(min_n = 2) |>
  normalise(method = "global_median") |>
  match_ptm_to_protein("id")

dpu_diagnostics(x)                      # which pairing do the data support?
ctr <- make_contrasts(x)
res <- run_methods(x, ctr, methods = "limma")
head(res[res$analysis == "ptm_DPU", ])  # usage = site minus protein
usage_beta_check(x, ctr)                # calls sensitive to the unit protein coefficient
write_results_xlsx(x, res, "results.xlsx")
cat(methods_paragraph(x, "limma"))
```

To compare the three model forms (requires MSstatsPTM and msqrob2):

```r
res3 <- run_methods(x, ctr, methods = c("limma", "msstatsptm", "msqrob2ptm"))
compare_methods(res3, fdr = 0.05, min_log2fc = 0)
```

One vignette per input software shows the import of real exports
(`browseVignettes("DiffPTM")` after installing with vignettes); small example exports of each
software are in `inst/extdata/formats/`.

## Statistical defaults, in brief

* **Usage model.** When the PTM-enriched and total-proteome samples share their names, each site
  is paired with its protein sample by sample (`blocking = "auto"`: the *per-site* paired model,
  equivalent to limma on the within-sample site-minus-protein differences). The consensus model
  (`blocking = "random"`) gains degrees of freedom but assumes that sites and proteins vary alike;
  `dpu_diagnostics()` reports when that holds.
* **Normalisation.** Each assay is centred on the global median of all samples
  (`method = "global_median"`), which assumes that the bulk of the features does not change
  between conditions. Per-condition centring (`method = "median"`) keeps a global shift between
  conditions and assumes instead that all conditions were loaded and enriched equally.
  `normalisation_diagnostics(x)`, run before `normalise()`, says which one the data support.
* **Missing values.** Tests use the observed values; imputation is optional.
* **Unit protein coefficient.** Every usage test assumes that, without a usage change, a site
  follows its protein one to one; `usage_beta_check()` flags calls that depend on it.

## Repository layout

```
R/               package code: importers, data model, processing, model forms, reports, interface
tests/testthat/  unit tests (simulated experiments with known truth, every importer)
vignettes/       interface walkthrough and one vignette per input software
inst/extdata/    example dataset and small example exports of each software
inst/benchmark/  simulation and spike-in benchmarks of the accompanying article
inst/app/        help page of the interface
Dockerfile       image serving the interface without installing R
```

## Citation

If you use DiffPTM, please cite the package (`citation("DiffPTM")`); the article is in preparation.

## Contributing and support

Bug reports and feature requests: [issues](https://github.com/Jacky11/DiffPTM/issues).
Please include the output of `sessionInfo()` and, if possible, a small example that
reproduces the problem.

## License

GPL-3. See [LICENSE.md](LICENSE.md).
