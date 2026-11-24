#################################################
# Analysis parameters
# One named list describes an analysis completely. It is filled either by the
# Shiny interface (server.R) or by a YAML file (scripts/run_analysis.R), and it
# is saved as parameters.yml next to the results so that any analysis can be
# re-run identically.
#################################################

default_params <- function() {
  list(
    input_file          = NULL,       # proteinGroups.txt (MaxQuant) or Spectronaut report
    spectronaut         = FALSE,      # TRUE for a Spectronaut report
    intensity           = "Raw",      # "Raw" (Intensity) or "LFQ" (LFQ intensity)
    delete_columns      = NULL,       # columns flagging with "+" the proteins to delete;
                                      #   NULL = all such columns (default of the app)
    cond1_name          = "Condition 1",
    cond1_samples       = NULL,       # column names, e.g. "Intensity.C_1B"
    cond2_name          = "Condition 2",
    cond2_samples       = NULL,
    ibaq                = TRUE,       # average iBAQ of present/absent proteins
    min_unique_peptides = 1,          # a value is quantified if >= this number of unique peptides
    min_observed        = 2,          # min. number of quantified values for a protein to be sufficiently quantified in a condition
    min_percent_observed = 0,         # min. % of quantified values (of the samples of the condition), combined with min_observed
    quantified_in       = "one",      # tested if sufficiently quantified in "one" (at least one) or "each" condition
    present_min         = NA,         # present/absent: min. quantified values where present (NA = same as above)
    absent_max          = 0,          # present/absent: max. quantified values where absent
    normalization       = "Median",   # "Median" or "None"
    imputation          = "MLE",      # "MLE", "SLSA" or "None"
    fold_change         = 2,          # minimal fold-change (not log2)
    test                = "LIMMA",    # "LIMMA" or "paired.LIMMA"
    pi0_method          = "pounds",   # "pounds" or "slim" (cp4p)
    conditional_pvalues = FALSE,      # p-values conditional on |log2FC| > threshold before FDR adjustment
    fdr                 = 0.01,
    saint               = FALSE,      # run SAINTexpress (cond1 = controls, cond2 = baits)
    saint_avgp          = 0.5,
    saint_exe           = NULL,       # NULL = bin/SAINTexpress-int(.exe)
    editable_pptx       = FALSE,      # extra PPTX with editable plots (package 'export')
    overlap             = FALSE,      # sample overlap (Jaccard / overlap indices) + pvclust AU/BP
    overlap_prefix      = "AllProteins",
    overlap_groups      = "",         # comma-separated group of each sample ("" = sample name without trailing number)
    overlap_nboot       = 10000,
    overlap_nperm       = 9999,
    n_cores             = 5,          # cores for pvclust (pvclust results depend on this number)
    seed                = 1234        # random seed (SLSA imputation; MLE imputation and overlap analyses use their own fixed seeds)
  )
}

# Read a YAML parameter file; relative paths are relative to the YAML file
read_params <- function(file) {
  p <- yaml::read_yaml(file)
  p <- utils::modifyList(default_params(), p)
  base <- dirname(normalizePath(file))
  if (!is.null(p$input_file) && !file.exists(p$input_file))
    p$input_file <- normalizePath(file.path(base, p$input_file), mustWork = FALSE)
  check_params(p)
}

write_params <- function(p, file) {
  yaml::write_yaml(p[names(p) != "data"], file)
}

check_params <- function(p) {
  p <- utils::modifyList(default_params(), p)
  stop_if <- function(cond, ...) if (cond) stop(..., call. = FALSE)
  stop_if(!p$intensity %in% c("Raw", "LFQ"), "intensity must be 'Raw' or 'LFQ'")
  stop_if(!p$quantified_in %in% c("one", "each"), "quantified_in must be 'one' or 'each'")
  stop_if(!is.numeric(p$min_observed) || p$min_observed < 1, "min_observed must be >= 1")
  stop_if(!is.numeric(p$min_percent_observed) || p$min_percent_observed < 0 || p$min_percent_observed > 100,
          "min_percent_observed must be between 0 and 100")
  stop_if(!is.numeric(p$absent_max) || p$absent_max < 0, "absent_max must be >= 0")
  stop_if(!is.null(p$present_min) && !is.na(p$present_min) && p$present_min < 1, "present_min must be >= 1")
  stop_if(!p$normalization %in% c("Median", "None"), "normalization must be 'Median' or 'None'")
  stop_if(!p$imputation %in% c("MLE", "SLSA", "None"), "imputation must be 'MLE', 'SLSA' or 'None'")
  stop_if(!p$test %in% c("LIMMA", "paired.LIMMA"), "test must be 'LIMMA' or 'paired.LIMMA'")
  stop_if(!p$pi0_method %in% c("pounds", "slim"), "pi0_method must be 'pounds' or 'slim'")
  stop_if(length(p$cond1_samples) < 2 || length(p$cond2_samples) < 2,
          "At least 2 samples are needed in each condition.")
  stop_if(identical(p$cond1_name, p$cond2_name), "The two conditions must have different names.")
  stop_if(length(intersect(p$cond1_samples, p$cond2_samples)) > 0,
          "A sample cannot be in both conditions: ", paste(intersect(p$cond1_samples, p$cond2_samples), collapse = ", "))
  stop_if(p$test == "paired.LIMMA" && length(p$cond1_samples) != length(p$cond2_samples),
          "The paired test needs the same number of samples in both conditions.")
  stop_if(!is.numeric(p$fdr) || p$fdr <= 0 || p$fdr >= 1, "fdr must be in ]0, 1[")
  stop_if(!is.numeric(p$fold_change) || p$fold_change < 1, "fold_change must be >= 1")
  stop_if(isTRUE(p$conditional_pvalues) && p$fold_change <= 1,
          "conditional_pvalues needs a fold_change threshold > 1")
  stop_if(isTRUE(p$editable_pptx) && !requireNamespace("export", quietly = TRUE),
          "editable_pptx needs the R package 'export': install.packages('export')")
  p
}

# Table "ParametersAnalysis" of the Excel file
parameters_table <- function(p, delete_columns) {
  param_label = c('Kind of intensities used for the analysis:',
                  'Columns indicating by + proteins to delete:',
                  'Name of the first condition:',
                  'Choose variables containing intensities for the 1st condition:',
                  'Name of the second condition:',
                  'Choose variables containing intensities for the 2nd condition:',
                  'IBAQ values:',
                  'Apply a present/absent filter:',
                  'Minimal number of unique identified peptides for a quantified value:',
                  'Minimal number of quantified values for a protein to be sufficiently quantified in a condition:',
                  'Minimal percentage of quantified values for a protein to be sufficiently quantified in a condition (between 0 and 100):',
                  'Normalization method:',
                  'Imputation method:',
                  'Select a minimal ratio (fold-change) between conditions (2 by default, i.e. log2(fold-change) threshold is fixed at log2(2)=1):',
                  'Statistical test:',
                  'Estimator of the proportion of differentially abundant proteins:',
                  'FDR level:',
                  'SAINT',
                  'Conditional p-values (given |log2FC| > threshold):',
                  'Proteins tested if sufficiently quantified in:',
                  'Present/absent: minimal number of quantified values in the condition where the protein is present:',
                  'Present/absent: maximal number of quantified values in the condition where the protein is absent:')
  param_values = c(p$intensity,
                   paste(delete_columns, collapse = "; "),
                   p$cond1_name,
                   paste(p$cond1_samples, collapse = "; "),
                   p$cond2_name,
                   paste(p$cond2_samples, collapse = "; "),
                   if (isTRUE(p$ibaq)) "Yes" else "No",
                   "YesPA",
                   p$min_unique_peptides,
                   p$min_observed,
                   p$min_percent_observed,
                   p$normalization,
                   p$imputation,
                   p$fold_change,
                   p$test,
                   p$pi0_method,
                   p$fdr,
                   p$saint,
                   isTRUE(p$conditional_pvalues),
                   if (p$quantified_in == "each") "each condition" else "at least one condition",
                   if (is.null(p$present_min) || is.na(p$present_min)) "same as above" else p$present_min,
                   p$absent_max)
  data.frame(cbind(param_label, param_values))
}
