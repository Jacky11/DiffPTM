#' Robustness of usage calls to the unit protein coefficient
#'
#' Every usage test of DiffPTM (and of MSstatsPTM and msqrob2PTM) assumes that, without a usage
#' change, a site follows its protein one to one (`site change = beta x protein change` with
#' `beta = 1`). When the true coefficient differs, the correction attributes `(beta - 1)` times
#' the protein change to usage, and protein-driven sites with a large protein change become false
#' usage calls. `beta` cannot be estimated per site from a few samples per condition, but the
#' sensitivity of each call to it can be checked: the per-site paired single model (or any model
#' of [method_limma()] selected through `...`) is refitted on `site - beta * protein` for each
#' value of `beta`, and a call is reported as robust when it stays significant for all of them.
#' The protein change of each site is recovered from the two extreme fits.
#'
#' @param x a `ptm_experiment` after [match_ptm_to_protein()].
#' @param contrasts contrasts, see [make_contrasts()].
#' @param beta coefficients to check (default 0.8 and 1.2, a 20 % departure).
#' @param fdr adjusted p-value threshold of a call.
#' @param ... passed to [method_limma()] (e.g. `blocking`, `pi0_method`).
#' @return a data frame with one row per site and contrast: usage log2FC and adjusted p-value
#'   at `beta = 1`, the protein log2 fold change, whether the site is called at each `beta`,
#'   and `robust_to_beta` (called at `beta = 1` and at every value of `beta`).
#' @examples
#' sim <- simulate_ptm_features(n_rep = 3, n_proteins = 300, n_sites = 200, seed = 1)
#' x <- match_ptm_to_protein(filter_min_observed(sim_to_experiment(sim), min_n = 2), "id")
#' head(usage_beta_check(x, make_contrasts(x), pi0_method = "bh"))
#' @export
usage_beta_check <- function(x, contrasts, beta = c(0.8, 1.2), fdr = 0.05, ...) {
  if (is.null(x$ptm$protein_intensity)) stop("Run match_ptm_to_protein() first")
  fit_b <- function(b) {
    xb <- x; xb$ptm$protein_intensity <- b * x$ptm$protein_intensity
    r <- method_limma(xb, contrasts, analyses = "ptm_DPU", ...)
    r[, c("feature_id", "contrast", "log2FC", "adj_pvalue")]
  }
  base <- fit_b(1)
  out <- data.frame(feature_id = base$feature_id, contrast = base$contrast, usage_log2FC = base$log2FC,
                    adj_pvalue = base$adj_pvalue, called = !is.na(base$adj_pvalue) & base$adj_pvalue < fdr)
  key <- paste(base$feature_id, base$contrast)
  fc <- list()
  for (b in beta) {
    r <- fit_b(b); m <- match(key, paste(r$feature_id, r$contrast))
    out[[paste0("called_beta_", b)]] <- !is.na(r$adj_pvalue[m]) & r$adj_pvalue[m] < fdr
    fc[[as.character(b)]] <- r$log2FC[m]
  }
  if (length(beta) >= 2) {                         # usage(beta) = site change - beta * protein change
    b1 <- min(beta); b2 <- max(beta)
    out$protein_log2FC <- (fc[[as.character(b1)]] - fc[[as.character(b2)]]) / (b2 - b1)
  }
  out$robust_to_beta <- out$called & Reduce(`&`, out[paste0("called_beta_", beta)])
  out
}
