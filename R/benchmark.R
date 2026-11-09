#' Simulate a PTM + proteome experiment with known truth
#'
#' Generates a `ptm_experiment` for benchmarking the DPA / DPU methods. Each PTM feature
#' belongs to one protein; its log2 intensity is the protein log2 intensity plus an offset
#' plus a site-specific usage effect, so the four classical situations are covered with
#' known labels:
#'
#' | class | protein changes | PTM changes | DPA truth | DPU truth |
#' |---|---|---|---|---|
#' | `null` | no | no | 0 | 0 |
#' | `protein_driven` | yes | follows protein | 1 | 0 |
#' | `usage_driven` | no | yes | 1 | 1 |
#' | `compensated` | yes | opposite usage change masks it | 0 | 1 |
#'
#' Within-sample PTM/protein correlation is controlled by `rho` (fraction of the sample
#' level variance shared by the PTM and its protein). Missing values are generated as a
#' mixture of MCAR and intensity-dependent MNAR (logistic in the intensity), separately
#' for the two assays.
#'
#' @param n_proteins number of proteins in the proteome.
#' @param n_ptm number of PTM features (parents drawn with replacement among proteins).
#' @param n_rep replicates per condition (two conditions, A and B).
#' @param effect log2 fold change applied to altered features.
#' @param frac fractions of the PTM features in classes `protein_driven`, `usage_driven`,
#'   `compensated` (the rest is `null`).
#' @param sd_bio biological standard deviation (log2) of protein and PTM measurements.
#' @param rho within-sample correlation between a PTM and its protein (0-1).
#' @param missing overall missing-value fraction per assay (`c(protein, ptm)`).
#' @param mnar fraction of missing values that are intensity-dependent (MNAR).
#' @param seed random seed.
#' @return a `ptm_experiment` with a `truth` attribute: data.frame(feature_id, class,
#'   dpa, dpu, true_log2fc_dpa, true_log2fc_dpu) and `protein_truth` (feature_id, changed,
#'   true_log2fc).
#' @export
simulate_ptm_experiment <- function(n_proteins = 2000, n_ptm = 1500, n_rep = 3, effect = 1,
                                    frac = c(protein_driven = 0.1, usage_driven = 0.1, compensated = 0.05),
                                    sd_bio = 0.3, rho = 0.5, missing = c(0.05, 0.2), mnar = 0.7, seed = 1) {
  set.seed(seed)
  cond <- rep(c("A", "B"), each = n_rep); samples <- paste0(cond, "_", seq_len(n_rep))
  n_s <- length(samples)
  # proteins
  parent <- sample(n_proteins, n_ptm, replace = TRUE)
  cls <- sample(c(rep("protein_driven", round(frac[["protein_driven"]] * n_ptm)),
                  rep("usage_driven",   round(frac[["usage_driven"]] * n_ptm)),
                  rep("compensated",    round(frac[["compensated"]] * n_ptm))))
  cls <- c(cls, rep("null", n_ptm - length(cls)))[seq_len(n_ptm)]
  prot_changed <- logical(n_proteins); prot_changed[parent[cls %in% c("protein_driven", "compensated")]] <- TRUE
  prot_changed[sample(setdiff(seq_len(n_proteins), parent), round(0.05 * n_proteins))] <- TRUE  # some changed proteins without PTM
  prot_fc <- ifelse(prot_changed, sample(c(-1, 1), n_proteins, TRUE) * effect, 0)
  base_p <- stats::rnorm(n_proteins, 24, 2.5)
  shared <- matrix(stats::rnorm(n_proteins * n_s, 0, sd_bio * sqrt(rho)), n_proteins)      # part shared with PTM
  own_p  <- matrix(stats::rnorm(n_proteins * n_s, 0, sd_bio * sqrt(1 - rho)), n_proteins)
  P <- base_p + shared + own_p + outer(prot_fc, as.numeric(cond == "B"))
  dimnames(P) <- list(paste0("PROT", seq_len(n_proteins)), samples)
  # PTMs
  usage_fc <- ifelse(cls == "usage_driven", sample(c(-1, 1), n_ptm, TRUE) * effect,
              ifelse(cls == "compensated", -prot_fc[parent], 0))
  own_q <- matrix(stats::rnorm(n_ptm * n_s, 0, sd_bio * sqrt(1 - rho)), n_ptm)
  Q <- base_p[parent] - 3 + shared[parent, , drop = FALSE] + own_q +
       outer(prot_fc[parent] + usage_fc, as.numeric(cond == "B"))
  dimnames(Q) <- list(paste0("PROT", parent, "_S", seq_len(n_ptm)), samples)
  # missing values: MCAR + MNAR (logistic in intensity, calibrated to reach the target rate)
  add_missing <- function(M, rate) {
    if (rate <= 0) return(M)
    n_na <- round(rate * length(M))
    n_mnar <- round(mnar * n_na); n_mcar <- n_na - n_mnar
    w <- stats::plogis(-(M - stats::quantile(M, 0.2)) / 1)           # low intensities more likely missing
    idx_mnar <- sample(length(M), n_mnar, prob = w)
    idx_mcar <- sample(setdiff(seq_along(M), idx_mnar), n_mcar)
    M[c(idx_mnar, idx_mcar)] <- NA
    M
  }
  P <- add_missing(P, missing[1]); Q <- add_missing(Q, missing[2])
  prot <- new_assay(P, data.frame(feature_id = rownames(P), protein_ids = rownames(P), gene = NA_character_))
  ptm  <- new_assay(Q, data.frame(feature_id = rownames(Q), protein_ids = paste0("PROT", parent),
                                  stripped_sequence = NA_character_, modified_sequence = NA_character_,
                                  residue = "S", position = as.character(seq_len(n_ptm)),
                                  localization_prob = 1, ptm_type = "Phospho"))
  ann <- make_annotation(samples, cond)
  x <- new_ptm_experiment(prot, ptm, annotation = list(protein = ann, ptm = ann),
                          params = list(simulation = list(n_rep = n_rep, effect = effect, rho = rho, missing = missing, mnar = mnar, seed = seed)))
  attr(x, "truth") <- data.frame(feature_id = rownames(Q), class = cls,
                                 dpa = as.integer((prot_fc[parent] + usage_fc) != 0),
                                 dpu = as.integer(usage_fc != 0),
                                 true_log2fc_dpa = -(prot_fc[parent] + usage_fc),   # contrast is A - B
                                 true_log2fc_dpu = -usage_fc, stringsAsFactors = FALSE)
  attr(x, "protein_truth") <- data.frame(feature_id = rownames(P), changed = as.integer(prot_changed),
                                         true_log2fc = -prot_fc, stringsAsFactors = FALSE)
  x
}

#' Evaluate method results against simulated truth
#'
#' Builds the truth table of a simulated experiment ([simulate_ptm_experiment()] or
#' [sim_to_experiment()]) and calls [evaluate_against_truth()].
#'
#' @param res tidy results from [run_methods()] on a simulated experiment.
#' @param x the simulated `ptm_experiment` (carries the `truth` attributes).
#' @param fdr nominal FDR threshold(s) at which calls are evaluated.
#' @param ... passed to [evaluate_against_truth()].
#' @return see [evaluate_against_truth()].
#' @export
evaluate_results <- function(res, x, fdr = c(0.01, 0.05), ...) {
  evaluate_against_truth(res, truth_table(x), fdr, ...)
}

#' Truth table of a simulated experiment
#' @param x simulated `ptm_experiment`.
#' @return data.frame(analysis, feature_id, positive, true_log2fc, class).
#' @export
truth_table <- function(x) {
  truth <- attr(x, "truth"); ptruth <- attr(x, "protein_truth")
  if (is.null(truth)) stop("`x` carries no simulation truth")
  rbind(
    data.frame(analysis = "protein_DPA", feature_id = ptruth$feature_id, positive = ptruth$changed,
               true_log2fc = ptruth$true_log2fc, class = ifelse(ptruth$changed == 1, "changed", "null")),
    data.frame(analysis = "ptm_DPA", feature_id = truth$feature_id, positive = truth$dpa,
               true_log2fc = truth$true_log2fc_dpa, class = truth$class),
    data.frame(analysis = "ptm_DPU", feature_id = truth$feature_id, positive = truth$dpu,
               true_log2fc = truth$true_log2fc_dpu, class = truth$class))
}

# one method x analysis (x contrast) block, merged on the truth side: untested truth rows are
# kept with NA statistics
merge_truth <- function(res, truth, keys) {
  out <- list()
  for (e in unique(res$method)) {
    r <- res[res$method == e, ]
    tr <- truth[truth$analysis %in% unique(r$analysis), ]
    m <- merge(r[, setdiff(names(r), "method")], tr, by = keys, all.y = TRUE)
    m$method <- e
    out[[e]] <- m
  }
  do.call(rbind, out)
}

roc_metrics <- function(p, positive) {
  p[is.na(p)] <- 1
  n1 <- sum(positive == 1); n0 <- sum(positive == 0)
  if (!n1 || !n0) return(c(auc = NA_real_, pauc = NA_real_, pr_auc = NA_real_))
  rk <- rank(-p)
  auc <- (sum(rk[positive == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  o <- order(p)                                              # ties keep their input order (neutral)
  tp <- cumsum(positive[o] == 1); fp <- cumsum(positive[o] == 0)
  tpr <- c(0, tp / n1); fpr <- c(0, fp / n0)
  keep <- fpr <= 0.05
  f <- c(fpr[keep], 0.05); t <- c(tpr[keep], stats::approx(fpr, tpr, 0.05, ties = max)$y)
  pauc <- sum(diff(f) * (utils::head(t, -1) + utils::tail(t, -1)) / 2) / 0.05
  prec <- tp / seq_along(tp)
  pr_auc <- sum(prec[positive[o] == 1]) / n1               # average precision
  c(auc = auc, pauc = pauc, pr_auc = pr_auc)
}

#' Evaluate method results against a truth table
#'
#' Every feature of the truth table counts: a feature an method did not test (filtered,
#' too many missing values, no parent protein) is "not called" and ranked last, so methods
#' that test fewer features are not favoured. Set `universe = "tested"` for the
#' tested-only definition used in earlier versions.
#'
#' @param res tidy results from [run_methods()].
#' @param truth data.frame with columns `analysis`, `feature_id`, `positive` (0/1),
#'   `true_log2fc`, optionally `class` and `contrast` (when the truth depends on the
#'   contrast, e.g. spike-in designs).
#' @param fdr nominal FDR threshold(s).
#' @param by_contrast also report metrics per contrast.
#' @param universe `"truth"` (default) or `"tested"`.
#' @return data.frame, one row per method x analysis (x contrast) x fdr:
#'   * counts `n_features`, `n_tested`, `n_called`, `TP`, `FP`, `FN`, `TN`;
#'   * at the threshold: `observed_fdr` (FDP), `FNR`, `FPR`, `sensitivity` (TPR),
#'     `specificity` (TNR), `precision`, `accuracy`, `F1`, `MCC`;
#'   * threshold-free: `auc` (ROC), `pauc` (ROC area for FPR <= 0.05, divided by 0.05),
#'     `pr_auc` (average precision);
#'   * calibration: `null_ks` (Kolmogorov-Smirnov distance of the tested null p-values to
#'     the uniform), `ci_coverage` (fraction of tested features whose true log2FC lies in
#'     log2FC +/- 1.96 SE);
#'   * effects: `fc_bias` and `fc_rmse` on tested true positives, `fc_rmse_null` on tested
#'     nulls.
#' @export
evaluate_against_truth <- function(res, truth, fdr = c(0.01, 0.05), by_contrast = FALSE,
                                   universe = c("truth", "tested")) {
  universe <- match.arg(universe)
  per_contrast <- "contrast" %in% names(truth)
  keys <- c("analysis", "feature_id", if (per_contrast) "contrast")
  m <- merge_truth(res, truth, keys)
  if (universe == "tested") m <- m[!is.na(m$adj_pvalue), ]
  metrics <- function(r, f) {
    tested <- !is.na(r$adj_pvalue); pos <- r$positive == 1
    called <- tested & r$adj_pvalue < f
    TP <- sum(called & pos); FP <- sum(called & !pos); FN <- sum(!called & pos); TN <- sum(!called & !pos)
    div <- function(a, b) if (b > 0) a / b else NA_real_
    mcc_den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    rm <- roc_metrics(r$pvalue, r$positive)
    se_ok <- tested & is.finite(r$se)
    nullp <- r$pvalue[tested & !pos]
    data.frame(fdr = f, n_features = nrow(r), n_tested = sum(tested), n_called = sum(called),
               TP = TP, FP = FP, FN = FN, TN = TN,
               observed_fdr = div(FP, TP + FP), FNR = div(FN, TP + FN), FPR = div(FP, FP + TN),
               sensitivity = div(TP, TP + FN), specificity = div(TN, TN + FP), precision = div(TP, TP + FP),
               accuracy = div(TP + TN, nrow(r)), F1 = div(2 * TP, 2 * TP + FP + FN),
               MCC = if (mcc_den > 0) (TP * TN - FP * FN) / mcc_den else NA_real_,
               auc = rm[["auc"]], pauc = rm[["pauc"]], pr_auc = rm[["pr_auc"]],
               null_ks = if (length(nullp) > 10) suppressWarnings(stats::ks.test(nullp, "punif")$statistic[[1]]) else NA_real_,
               ci_coverage = if (sum(se_ok)) mean(abs(r$log2FC[se_ok] - r$true_log2fc[se_ok]) <= stats::qnorm(0.975) * r$se[se_ok]) else NA_real_,
               fc_bias = if (sum(tested & pos)) mean((r$log2FC - r$true_log2fc)[tested & pos]) else NA_real_,
               fc_rmse = if (sum(tested & pos)) sqrt(mean((r$log2FC - r$true_log2fc)[tested & pos]^2)) else NA_real_,
               fc_rmse_null = if (sum(tested & !pos)) sqrt(mean((r$log2FC - r$true_log2fc)[tested & !pos]^2)) else NA_real_)
  }
  groups <- c("method", "analysis", if (by_contrast) "contrast")
  rows <- lapply(split(m, interaction(m[groups], drop = TRUE)), function(r)
    do.call(rbind, lapply(fdr, function(f) cbind(r[1, groups, drop = FALSE], metrics(r, f)))))
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}

#' Call rate per truth class
#'
#' For each method x analysis x nominal FDR and each truth class (`null`,
#' `protein_driven`, `usage_driven`, `compensated`), the fraction of features called.
#' For DPU, `protein_driven` sites are nulls: their call rate is the error that protein
#' correction is meant to remove; `compensated` sites are positives that abundance tests
#' cannot see.
#' @inheritParams evaluate_against_truth
#' @export
evaluate_by_class <- function(res, truth, fdr = c(0.01, 0.05)) {
  if (!"class" %in% names(truth)) stop("`truth` needs a `class` column")
  m <- merge_truth(res, truth, c("analysis", "feature_id"))
  rows <- list()
  for (f in fdr) {
    m$called <- !is.na(m$adj_pvalue) & m$adj_pvalue < f
    a <- stats::aggregate(cbind(call_rate = called, n = 1) ~ method + analysis + class, m,
                          function(v) if (length(v)) mean(v) else NA_real_)
    a$n <- as.integer(stats::aggregate(called ~ method + analysis + class, m, length)$called)
    a$fdr <- f
    rows[[length(rows) + 1]] <- a
  }
  do.call(rbind, rows)
}

#' Observed false discovery proportion over a range of nominal FDR levels
#' @inheritParams evaluate_against_truth
#' @param grid nominal FDR levels.
#' @return data.frame(method, analysis, nominal, observed_fdr, n_called).
#' @export
fdp_curve <- function(res, truth, grid = c(0.001, 0.005, seq(0.01, 0.2, by = 0.01))) {
  m <- merge_truth(res, truth, c("analysis", "feature_id"))
  out <- list()
  for (k in split(m, interaction(m$method, m$analysis, drop = TRUE))) {
    called_at <- function(f) !is.na(k$adj_pvalue) & k$adj_pvalue < f
    out[[length(out) + 1]] <- data.frame(method = k$method[1], analysis = k$analysis[1], nominal = grid,
      observed_fdr = vapply(grid, function(f) { c <- called_at(f); if (sum(c)) mean(k$positive[c] == 0) else NA_real_ }, 0),
      n_called = vapply(grid, function(f) sum(called_at(f)), 0))
  }
  do.call(rbind, out)
}
