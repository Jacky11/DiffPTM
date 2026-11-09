#' limma method: DPA on proteins and PTM features, DPU with the single-model approach
#'
#' * `protein_DPA` / `ptm_DPA`: `lmFit` on the assay with a no-intercept group design,
#'   contrasts `group1 - group2`, `eBayes(robust = TRUE)`.
#' * `ptm_DPU`: the PTM intensities and the matched parent-protein intensities are stacked
#'   column-wise for each PTM feature; the group factor is `component:group`
#'   (component = ptm / prot) and the tested contrast is
#'   `(ptm.g1 - ptm.g2) - (prot.g1 - prot.g2)`, i.e. the interaction term of the
#'   single-model strategy (Wolski et al. 2027, §1.5). Missing values are tolerated by
#'   `lmFit`, so no complete-case filtering is applied; features with fewer than
#'   `min_obs_per_group` observed values in a tested group are reported with `NA`.
#'
#' Pairing / blocking: when the annotation has a `subject` column, `blocking = "fixed"`
#' adds the subject as a fixed effect (what the original paired tests did);
#' `"random"` uses `duplicateCorrelation`; `"auto"` picks `"fixed"` when every subject
#' appears in every group, `"random"` otherwise.
#'
#' For the DPU model, when the PTM and proteome assays were measured on the same samples
#' (identical sample names) and no `subject` is given, `pair_by_sample = TRUE` (default)
#' uses the sample as the blocking factor of the stacked model. This captures the
#' within-sample variation shared by a PTM and its parent protein, which the unpaired
#' single model and MSstatsPTM ignore; without it the DPU test is strongly conservative
#' whenever PTM and protein intensities are correlated within samples.
#'
#' Per-site pairing (default): with `blocking = "auto"` the sample enters the DPU model as a
#' fixed effect. Each sample term absorbs what its site and protein values have in common, so
#' the residuals are the within-sample differences: this is the exact paired analysis of
#' correct-then-model (as msqrob2PTM, without robust regression), with each site's own
#' covariance, its own residual variance of the difference, and `2n - 2` residual degrees of
#' freedom for two conditions of `n` samples. It assumes neither that the site and its protein
#' have the same residual variance nor that all sites share one correlation.
#'
#' Consensus pairing (`blocking = "random"`): the sample term is a random effect fitted with
#' `limma::duplicateCorrelation`, i.e. one consensus within-sample correlation estimated from
#' all sites; each site keeps the residual degrees of freedom of the unpaired model (`4n - 4`),
#' and its standard error is approximately the unpaired one times `sqrt(1 - rho_consensus)`.
#' This doubles the residual degrees of freedom, but only legitimately when the site and its
#' protein have similar residual variances and similar correlations across sites: the model
#' pools the two residual variances into one. In the real datasets of the DiffPTM benchmark
#' the residual variance of a site was typically 5-10 times that of its protein, varying
#' between sites, and in simulations with such variances the consensus model was
#' anti-conservative, whereas the per-site model stayed closest to the nominal level. Use
#' [dpu_diagnostics()] to check these conditions; the consensus model can be preferred with two
#' replicates per group when they hold.
#'
#' Complete pairs: with the consensus correlation, a sample in which only the site or only
#' the protein is quantified would still enter the fit, and the shared correlation would let
#' this unpaired value inform the other component; in the benchmark this made the
#' all-observation consensus fit (default of earlier development versions of DiffPTM) anti-conservative.
#' `complete_pairs = TRUE` (default) sets unpaired values to `NA` before fitting the consensus
#' model, so that only samples with both values enter. The per-site paired model ignores
#' unpaired values by construction.
#'
#' @inheritParams run_methods
#' @param ... ignored (lets [run_methods()] pass method-specific arguments).
#' @param robust passed to `limma::eBayes`.
#' @param blocking `"auto"`, `"none"`, `"fixed"`, `"random"`. For the DPU model with sample
#'   pairing, `"auto"` is the per-site paired model (`"fixed"`); for designs with a `subject`,
#'   `"auto"` picks `"fixed"` when every subject appears in every group and `"random"` otherwise.
#' @param min_obs_per_group minimal observed values per group for a contrast to be
#'   reported.
#' @param pair_by_sample see Details.
#' @param complete_pairs restrict the consensus-correlation DPU fit (`blocking = "random"`) to
#'   samples in which both the site and its protein are quantified (default; see Details). `FALSE` gives the
#'   all-observation consensus fit of earlier development versions of DiffPTM, which was anti-conservative for sites
#'   with unpaired values.
#' @param moderate use empirical Bayes variance moderation (`limma::eBayes`; default). With
#'   `FALSE`, ordinary t-statistics are computed from each feature's own residual variance
#'   (`limma::lmFit` only, residual degrees of freedom): used in the benchmark to isolate the
#'   effect of moderation from that of pairing.
#' @param weights use precision weights: the number of features (peptides, peptidoforms)
#'   behind each summarised value, when the assays carry it (`n_features`, e.g. from
#'   [summarise_features()]). A value summarised from more features is less noisy; without
#'   weights, limma assumes equal precision within a feature, which makes the paired DPU
#'   test anti-conservative when the site-protein correlation is high and the number of
#'   observed features varies between samples. Exploratory; default `FALSE`.
#' @export
method_limma <- function(x, contrasts, analyses = c("protein_DPA", "ptm_DPA", "ptm_DPU"),
                         pi0_method = "pounds", robust = TRUE,
                         blocking = c("auto", "none", "fixed", "random"), min_obs_per_group = 2,
                         pair_by_sample = TRUE, weights = FALSE, moderate = TRUE, complete_pairs = TRUE, ...) {
  blocking <- match.arg(blocking)
  res <- list()
  if ("protein_DPA" %in% analyses && !is.null(x$protein)) {
    ann <- x$annotation$protein
    W <- if (weights) x$protein$n_features
    res$protein_DPA <- fit_limma(x$protein$intensity, groups_of(x, "protein"), ann$subject, mixture = ann$mixture, weights = W,
                                 contrasts, blocking, robust, pi0_method, min_obs_per_group,
                                 analysis = "protein_DPA", moderate = moderate)
  }
  if ("ptm_DPA" %in% analyses && !is.null(x$ptm)) {
    ann <- x$annotation$ptm
    W <- if (weights) x$ptm$n_features
    res$ptm_DPA <- fit_limma(x$ptm$intensity, groups_of(x, "ptm"), ann$subject, mixture = ann$mixture, weights = W,
                             contrasts, blocking, robust, pi0_method, min_obs_per_group,
                             analysis = "ptm_DPA", moderate = moderate)
  }
  if ("ptm_DPU" %in% analyses && !is.null(x$ptm)) {
    if (is.null(x$ptm$protein_intensity)) stop("Run match_ptm_to_protein() before the DPU analysis")
    S_dpu <- x$ptm$intensity; P_dpu <- x$ptm$protein_intensity
    by_sample <- pair_by_sample && all(is.na(c(x$annotation$ptm$subject, x$annotation$protein$subject))) &&
      setequal(colnames(S_dpu), colnames(P_dpu)) && blocking == "random"
    if (complete_pairs && by_sample) {
      Pa <- P_dpu[, colnames(S_dpu), drop = FALSE]; unpaired <- is.na(S_dpu) | is.na(Pa)
      S_dpu[unpaired] <- NA; Pa[unpaired] <- NA; P_dpu[, colnames(S_dpu)] <- Pa
    }
    M <- cbind(S_dpu, P_dpu)
    comp <- c(rep("ptm", ncol(x$ptm$intensity)), rep("prot", ncol(x$ptm$protein_intensity)))
    grp <- c(as.character(groups_of(x, "ptm")), as.character(groups_of(x, "protein")))
    subj <- c(x$annotation$ptm$subject, x$annotation$protein$subject)
    same_samples <- setequal(colnames(x$ptm$intensity), colnames(x$ptm$protein_intensity))
    if (pair_by_sample && all(is.na(subj)) && same_samples && blocking != "none")
      subj <- c(colnames(x$ptm$intensity), colnames(x$ptm$protein_intensity))
    # sample pairing: "auto" is the per-site paired model (fixed sample term)
    dpu_blocking <- if (pair_by_sample && all(is.na(c(x$annotation$ptm$subject, x$annotation$protein$subject))) &&
                        same_samples && blocking == "auto") "fixed" else blocking
    stacked_group <- factor(paste(comp, grp, sep = "__"))
    mix <- c(x$annotation$ptm$mixture, x$annotation$protein$mixture)
    # each fraction has its own mixture (TMT plex) effect: the enriched and total fractions are labelled
    # and measured separately, so a plex offset of the sites need not be shared by the proteins
    if (!all(is.na(mix))) mix <- paste(comp, mix, sep = "_")
    W <- if (weights && !is.null(x$ptm$n_features) && !is.null(x$ptm$protein_n_features))
      cbind(x$ptm$n_features, x$ptm$protein_n_features)
    res$ptm_DPU <- fit_limma(M, stacked_group, subj, mixture = mix, weights = W, contrasts, dpu_blocking, robust, pi0_method,
                             min_obs_per_group, analysis = "ptm_DPU", interaction = TRUE, moderate = moderate)
  }
  out <- do.call(rbind, res)
  if (is.null(out)) return(empty_result())
  out$method <- "limma"
  rownames(out) <- NULL
  out
}

#' @keywords internal
fit_limma <- function(M, group, subject, mixture = NULL, weights = NULL, contrasts, blocking, robust, pi0_method,
                      min_obs_per_group, analysis, interaction = FALSE, moderate = TRUE) {
  group <- factor(group)
  use_subject <- !is.null(subject) && !all(is.na(subject)) && blocking != "none"
  if (use_subject && blocking == "auto") {
    tab <- table(subject, group)
    blocking <- if (all(tab > 0)) "fixed" else "random"
  }
  design <- make_design(group, if (use_subject && blocking == "fixed") subject else NULL, mixture)
  lv_map <- attr(design, "group_levels")           # original level -> design column name

  # "A" -> design column of A; "B+C" -> (B + C)/2 ; `prefix` selects the stacked component
  gexpr <- function(s, prefix = "") {
    m <- paste0(prefix, trimws(strsplit(s, "+", fixed = TRUE)[[1]]))
    if (any(is.na(lv_map[m]))) stop("Unknown group in contrast: ", s)
    if (length(m) == 1) unname(lv_map[m]) else paste0("(", paste(lv_map[m], collapse = " + "), ")/", length(m))
  }
  exprs <- if (interaction) {
    mapply(function(g1, g2) paste0("(", gexpr(g1, "ptm__"), " - ", gexpr(g2, "ptm__"), ") - (",
                                   gexpr(g1, "prot__"), " - ", gexpr(g2, "prot__"), ")"),
           contrasts$group1, contrasts$group2, USE.NAMES = FALSE)
  } else {
    mapply(function(g1, g2) paste0(gexpr(g1), " - ", gexpr(g2)), contrasts$group1, contrasts$group2, USE.NAMES = FALSE)
  }
  cm <- limma::makeContrasts(contrasts = exprs, levels = design)
  colnames(cm) <- contrasts$name

  quiet_fit <- function(expr) withCallingHandlers(expr, warning = function(w) {
    if (grepl("Partial NA coefficients", conditionMessage(w))) invokeRestart("muffleWarning") })
  random_block <- use_subject && blocking == "random"
  corr <- if (random_block) limma::duplicateCorrelation(M, design, block = factor(subject))$consensus.correlation
  if (random_block && !is.finite(corr))
    stop("The consensus within-sample correlation could not be estimated (too few residual degrees of freedom per ",
         "site for this design); use the per-site paired model (blocking = \"fixed\" or \"auto\").", call. = FALSE)
  if (!is.null(weights)) { weights <- as.matrix(weights); weights[!is.finite(weights) | weights <= 0] <- 1 }
  do_fit <- function(X) quiet_fit(if (random_block) limma::lmFit(M, X, block = factor(subject), correlation = corr, weights = weights)
                                  else limma::lmFit(M, X, weights = weights))
  # limma::contrasts.fit() derives standard errors from the full design's covariance, which is
  # exact only for orthogonal designs or complete data. With blocking columns (sample,
  # subject, plex) and missing values it is not, and the DPU test becomes anti-conservative.
  # In that case each contrast is refitted as the first coefficient of a reparametrised
  # design (limma::contrastAsCoef), which gives exact per-feature standard errors.
  # moderated (eBayes) or ordinary statistics from a fitted MArrayLM object
  stats_of <- function(fit) {
    if (moderate) {
      fit <- limma::eBayes(fit, robust = robust)
      list(est = fit$coefficients, t = fit$t, p = fit$p.value, se = sqrt(fit$s2.post) * fit$stdev.unscaled,
           df = matrix(fit$df.total, nrow(fit$coefficients), ncol(fit$coefficients)))
    } else {
      se <- fit$sigma * fit$stdev.unscaled
      t <- fit$coefficients / se
      list(est = fit$coefficients, t = t, p = 2 * stats::pt(-abs(t), fit$df.residual), se = se,
           df = matrix(fit$df.residual, nrow(fit$coefficients), ncol(fit$coefficients)))
    }
  }
  exact <- (anyNA(M) || !is.null(weights)) && ncol(design) > nlevels(group)
  if (!exact) {
    st <- stats_of(limma::contrasts.fit(do_fit(design), cm))
    est <- st$est; t <- st$t; pv <- st$p; se <- st$se; dfm <- st$df
  } else {
    est <- t <- pv <- se <- dfm <- matrix(NA_real_, nrow(M), ncol(cm), dimnames = list(rownames(M), colnames(cm)))
    for (k in seq_len(ncol(cm))) {
      X <- limma::contrastAsCoef(design, cm[, k, drop = FALSE], first = TRUE)$design
      st <- stats_of(do_fit(X))
      est[, k] <- st$est[, 1]; t[, k] <- st$t[, 1]; pv[, k] <- st$p[, 1]; se[, k] <- st$se[, 1]; dfm[, k] <- st$df[, 1]
    }
  }

  obs <- !is.na(M)
  n_by_group <- sapply(levels(group), function(lv) rowSums(obs[, group == lv, drop = FALSE]))
  if (is.null(dim(n_by_group))) n_by_group <- matrix(n_by_group, ncol = nlevels(group), dimnames = list(NULL, levels(group)))

  out <- lapply(seq_len(ncol(cm)), function(k) {
    involved <- rownames(cm)[cm[, k] != 0]
    orig <- names(lv_map)[match(involved, lv_map)]
    ok <- rowSums(n_by_group[, orig, drop = FALSE] >= min_obs_per_group) == length(orig)
    p <- pv[, k]; p[!ok] <- NA
    adj <- adjust_pvalues(p, pi0_method)
    data.frame(feature_id = rownames(M), analysis = analysis, contrast = colnames(cm)[k],
               log2FC = ifelse(ok, est[, k], NA_real_),
               se = ifelse(ok, se[, k], NA_real_),
               df = ifelse(ok, dfm[, k], NA_real_),
               statistic = ifelse(ok, t[, k], NA_real_),
               pvalue = p, adj_pvalue = adj$adjusted, method = "limma", stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  attr(out, "pi0") <- NULL
  out
}
