#' Add missing values with the imp4p mechanism
#'
#' Reproduces the missing-value process of `imp4p::sim.data()` on any complete log2
#' matrix (columns = samples). In each column, `round(rate * nrow)` values are removed:
#' a fraction `pi_mcar` uniformly at random (MCAR), the rest (MNAR) by weighted sampling
#' without replacement among the remaining values with probability
#' `max(1 - para * (x - xmin) / (xmax - xmin), 1e-8)`, where `xmin` / `xmax` are the 1 %
#' and 99 % quantiles of the column (imp4p uses the normal quantiles of the column's mean
#' structure; the empirical quantiles are their direct analogue for arbitrary data).
#'
#' @param M complete numeric matrix (log2 scale).
#' @param rate fraction of missing values per column.
#' @param pi_mcar fraction of the missing values that are MCAR.
#' @param para MNAR shape parameter (imp4p default 3; larger = MNAR concentrated at low
#'   intensities).
#' @return `M` with `NA`s; attribute `mcar` is a logical matrix flagging MCAR cells.
#' @export
add_missing_imp4p <- function(M, rate = 0.2, pi_mcar = 0.2, para = 3) {
  mcar <- matrix(FALSE, nrow(M), ncol(M))
  n_miss <- round(rate * nrow(M)); n_mcar <- floor(n_miss * pi_mcar); n_mnar <- n_miss - n_mcar
  for (j in seq_len(ncol(M))) {
    x <- M[, j]
    if (n_mcar > 0) {
      i_mcar <- sample.int(nrow(M), n_mcar)
      x[i_mcar] <- NA; mcar[i_mcar, j] <- TRUE
    }
    if (n_mnar > 0) {
      q <- stats::quantile(M[, j], c(0.01, 0.99), names = FALSE)
      p <- pmax(1 - para * (M[, j] - q[1]) / (q[2] - q[1]), 1e-8)
      avail <- which(!is.na(x))
      i_mnar <- avail[sample.int(length(avail), min(n_mnar, length(avail)), prob = p[avail])]
      x[i_mnar] <- NA
    }
    M[, j] <- x
  }
  attr(M, "mcar") <- mcar
  M
}

#' Simulate a feature-level PTM + proteome experiment with known truth
#'
#' Peptide-level total proteome and peptidoform-level PTM data, following the
#' hierarchical variance structure of `imp4p::sim.data()` (protein means `N(m_c, sd_c)`,
#' biological-sample effects `N(0, sd_rb)`, measurement noise `N(0, sd_r)`), with
#' differential effects added between two conditions A and B.
#'
#' * Proteins have `1 + NB(mu = 4, size = 2)` peptides (capped at 20), each with an
#'   ionisation offset `N(0, 1)`.
#' * Each site belongs to a random protein; it is measured by 1-3 peptidoforms. Its
#'   biological-sample effect is `rho * u + sqrt(1 - rho^2) * w`, where `u` is the
#'   parent protein's sample effect and `w` an independent `N(0, sd_rb)`: `rho` is the
#'   within-sample correlation between a site and its protein.
#' * Truth classes as in [simulate_ptm_experiment()]: `null`, `protein_driven` (DPA+ DPU-),
#'   `usage_driven` (DPA+ DPU+), `compensated` (DPA- DPU+); 5 % of the proteins without
#'   sites change too.
#' * Missing values: [add_missing_imp4p()] on each feature matrix.
#'
#' @param n_proteins,n_sites numbers of proteins and sites.
#' @param n_rep replicates per condition.
#' @param effect absolute log2 fold change of altered proteins / sites: one value, or a
#'   range `c(min, max)` from which each effect is drawn uniformly (sign at random).
#' @param frac fractions of sites in classes `protein_driven`, `usage_driven`, `compensated`.
#' @param m_c,sd_c mean and SD of protein log2 intensities (imp4p `m.c`, `sd.c`).
#' @param sd_rb biological-sample SD (imp4p `sd.rb`; default 0.3 instead of imp4p's 0.5,
#'   typical of cell-line biological replicates).
#' @param sd_r measurement SD of each feature (imp4p `sd.r`).
#' @param rho within-sample site-protein correlation: one value for all sites, or one value
#'   per site (vector of length `n_sites`) to simulate heterogeneous correlations.
#' @param beta coefficient relating a site's intensity to its protein's abundance change
#'   (1 by default: a site follows its protein one to one). A protein-driven site then changes
#'   by `beta` times the protein change, so that a unit-coefficient correction estimates a
#'   spurious usage change of `(beta - 1)` times the protein change; used in the sensitivity
#'   analysis of the unit-coefficient assumption.
#' @param rate,pi_mcar,para missing-value parameters, see [add_missing_imp4p()];
#'   `rate` may be a vector `c(protein, ptm)`.
#' @param seed random seed.
#' @param site_clustering 0 (default) assigns sites to proteins uniformly at random; a value
#'   `a > 0` draws protein weights from `Gamma(a)`, so that sites cluster on fewer proteins
#'   (`a = 0.15` gives a median of about three sites per site-carrying protein).
#' @param rho_condition optional `c(rho_A, rho_B)`: a within-sample correlation specific to
#'   each condition (overrides `rho`).
#' @param sd_heterogeneity SD of a log-normal factor multiplying the sample-effect SD of each
#'   protein and of each site (0: equal variances).
#' @param site_sd_factor multiplies the SD of the site's biological sample effects, keeping their
#'   correlation `rho` with the protein's (1: same SD as the protein; in real data the residual
#'   variance of a site is typically several times that of its protein, see the vignette).
#' @param batch_sd `c(protein, ptm)`: SD of feature-specific effects of two processing batches
#'   balanced within each condition and not given to the models (0: none).
#' @param loading_sd SD of a per-sample offset of the PTM assay only (enrichment efficiency),
#'   to be removed by normalisation.
#' @param outlier_frac,outlier_size fraction of site x sample values shifted by a random sign
#'   times `U(outlier_size)` (applied to all peptidoforms of the site in that sample).
#' @param rate_B_factor multiplies the PTM missing rate in condition B (unbalanced missingness).
#'   `n_rep` may be `c(n_A, n_B)` for unbalanced designs.
#' @return list with `protein_features`, `ptm_features` (log2 matrices, features x
#'   samples, with `NA`s), `protein_map` (feature -> protein), `ptm_map` (feature ->
#'   site, parent protein), `annotation`, `truth` (site-level), `protein_truth`.
#' @export
simulate_ptm_features <- function(n_proteins = 2000, n_sites = 1500, n_rep = 3, effect = c(0.5, 2),
                                  frac = c(protein_driven = 0.1, usage_driven = 0.1, compensated = 0.05),
                                  m_c = 25, sd_c = 2, sd_rb = 0.3, sd_r = 0.2, rho = 0.5,
                                  rate = c(0.1, 0.2), pi_mcar = 0.2, para = 3, seed = 1, beta = 1,
                                  site_clustering = 0, rho_condition = NULL, sd_heterogeneity = 0, batch_sd = c(0, 0),
                                  loading_sd = 0, outlier_frac = 0, outlier_size = c(1.5, 3), rate_B_factor = 1,
                                  site_sd_factor = 1) {
  force(rho); force(beta); force(rho_condition)   # evaluate arguments before seeding (lazy evaluation)
  # every option below draws random numbers only when it is switched on, so that the default
  # simulation is unchanged
  set.seed(seed)
  if (length(rate) == 1) rate <- c(rate, rate)
  if (!length(rho) %in% c(1, n_sites)) stop("rho must have length 1 or n_sites")
  draw_effect <- function(n) sample(c(-1, 1), n, TRUE) *
    (if (length(effect) == 1) rep(effect, n) else stats::runif(n, effect[1], effect[2]))
  if (length(n_rep) == 1) n_rep <- rep(n_rep, 2)
  cond <- rep(c("A", "B"), times = n_rep); samples <- paste0(cond, "_", c(seq_len(n_rep[1]), seq_len(n_rep[2]))); n_s <- length(samples)
  inB <- as.numeric(cond == "B")
  prot_ids <- sprintf("PROT%04d", seq_len(n_proteins))

  # site classes and protein effects
  parent <- if (site_clustering > 0) sample(n_proteins, n_sites, replace = TRUE, prob = stats::rgamma(n_proteins, shape = site_clustering)) else
    sample(n_proteins, n_sites, replace = TRUE)
  k <- round(frac * n_sites)
  cls <- sample(c(rep("protein_driven", k[["protein_driven"]]), rep("usage_driven", k[["usage_driven"]]),
                  rep("compensated", k[["compensated"]]), rep("null", n_sites - sum(k))))
  prot_changed <- logical(n_proteins)
  prot_changed[parent[cls %in% c("protein_driven", "compensated")]] <- TRUE
  free <- setdiff(seq_len(n_proteins), parent)
  prot_changed[free[sample.int(length(free), min(length(free), round(0.05 * n_proteins)))]] <- TRUE
  prot_fc <- ifelse(prot_changed, draw_effect(n_proteins), 0)
  usage_fc <- ifelse(cls == "usage_driven", draw_effect(n_sites),
              ifelse(cls == "compensated", -prot_fc[parent], 0))

  mu_p <- stats::rnorm(n_proteins, m_c, sd_c)
  u <- matrix(stats::rnorm(n_proteins * n_s, 0, sd_rb), n_proteins)          # protein sample effects
  if (sd_heterogeneity > 0) u <- u * exp(stats::rnorm(n_proteins, 0, sd_heterogeneity))

  # protein peptides
  n_pep <- pmin(1 + stats::rnbinom(n_proteins, mu = 4, size = 2), 20)
  pep_prot <- rep(seq_len(n_proteins), n_pep)
  pep_off <- stats::rnorm(length(pep_prot), 0, 1)
  P <- mu_p[pep_prot] + pep_off + u[pep_prot, , drop = FALSE] + outer(prot_fc[pep_prot], inB) +
       matrix(stats::rnorm(length(pep_prot) * n_s, 0, sd_r), length(pep_prot))
  pep_id <- paste0(prot_ids[pep_prot], "_pep", stats::ave(pep_prot, pep_prot, FUN = seq_along))
  dimnames(P) <- list(pep_id, samples)

  # sites and peptidoforms
  site_id <- sprintf("%s_S%d", prot_ids[parent], seq_len(n_sites))
  w <- matrix(stats::rnorm(n_sites * n_s, 0, sd_rb), n_sites)
  if (sd_heterogeneity > 0) w <- w * exp(stats::rnorm(n_sites, 0, sd_heterogeneity))
  rho_m <- if (is.null(rho_condition)) rho else matrix(rep(ifelse(cond == "A", rho_condition[1], rho_condition[2]), each = n_sites), n_sites)
  v <- rho_m * u[parent, , drop = FALSE] + sqrt(1 - rho_m^2) * w
  if (site_sd_factor != 1) v <- v * site_sd_factor
  site_mu <- mu_p[parent] - 3 + stats::rnorm(n_sites, 0, 1)
  n_form <- 1 + stats::rbinom(n_sites, 2, 0.3)
  form_site <- rep(seq_len(n_sites), n_form)
  Q <- site_mu[form_site] + stats::rnorm(length(form_site), 0, 0.7) + v[form_site, , drop = FALSE] +
       outer(beta * prot_fc[parent[form_site]] + usage_fc[form_site], inB) +
       matrix(stats::rnorm(length(form_site) * n_s, 0, sd_r), length(form_site))
  form_id <- paste0(site_id[form_site], "_f", stats::ave(form_site, form_site, FUN = seq_along))
  dimnames(Q) <- list(form_id, samples)

  if (outlier_frac > 0) {
    n_out <- round(outlier_frac * n_sites * n_s); cell <- sample.int(n_sites * n_s, n_out)
    shift <- matrix(0, n_sites, n_s); shift[cell] <- sample(c(-1, 1), n_out, TRUE) * stats::runif(n_out, outlier_size[1], outlier_size[2])
    Q <- Q + shift[form_site, , drop = FALSE]
  }
  if (any(batch_sd > 0)) {
    batch <- stats::ave(seq_along(cond), cond, FUN = function(i) (seq_along(i) - 1) %% 2 + 1)
    P <- P + matrix(stats::rnorm(nrow(P) * 2, 0, batch_sd[1]), nrow(P))[, batch, drop = FALSE]
    Q <- Q + matrix(stats::rnorm(nrow(Q) * 2, 0, batch_sd[2]), nrow(Q))[, batch, drop = FALSE]
  }
  if (loading_sd > 0) Q <- sweep(Q, 2, stats::rnorm(n_s, 0, loading_sd), "+")
  P <- add_missing_imp4p(P, rate[1], pi_mcar, para)
  if (rate_B_factor == 1) Q <- add_missing_imp4p(Q, rate[2], pi_mcar, para) else {
    jA <- cond == "A"
    Q[, jA] <- add_missing_imp4p(Q[, jA, drop = FALSE], rate[2], pi_mcar, para)
    Q[, !jA] <- add_missing_imp4p(Q[, !jA, drop = FALSE], min(rate[2] * rate_B_factor, 0.9), pi_mcar, para)
  }

  # a site drawn as null or usage-driven can sit on a protein changed by another of its sites:
  # the class is therefore recomputed from the effects actually applied
  dpa <- as.integer((beta * prot_fc[parent] + usage_fc) != 0); dpu <- as.integer(usage_fc != 0)
  cls <- ifelse(dpu == 0, ifelse(dpa == 0, "null", "protein_driven"), ifelse(dpa == 1, "usage_driven", "compensated"))
  list(
    protein_features = P, ptm_features = Q,
    protein_map = data.frame(feature = pep_id, protein = prot_ids[pep_prot], stringsAsFactors = FALSE),
    ptm_map = data.frame(feature = form_id, site = site_id[form_site], protein = prot_ids[parent[form_site]],
                         stringsAsFactors = FALSE),
    annotation = make_annotation(samples, cond),
    truth = data.frame(feature_id = site_id, protein = prot_ids[parent], class = cls, dpa = dpa, dpu = dpu,
                       true_log2fc_dpa = -(beta * prot_fc[parent] + usage_fc), true_log2fc_dpu = -usage_fc,
                       stringsAsFactors = FALSE),
    protein_truth = data.frame(feature_id = prot_ids, changed = as.integer(prot_changed),
                               true_log2fc = -prot_fc, stringsAsFactors = FALSE),
    params = list(n_rep = n_rep, effect = effect, rho = rho, rate = rate, pi_mcar = pi_mcar, para = para,
                  sd_rb = sd_rb, sd_r = sd_r, seed = seed, beta = beta, site_clustering = site_clustering,
                  rho_condition = rho_condition, sd_heterogeneity = sd_heterogeneity, batch_sd = batch_sd,
                  loading_sd = loading_sd, outlier_frac = outlier_frac, rate_B_factor = rate_B_factor,
                  site_sd_factor = site_sd_factor))
}

#' Summarise features to proteins or sites by Tukey median polish
#'
#' Per group of features (rows of a log2 matrix), fits `stats::medpolish()` on the
#' observed values and returns overall + column effect for each sample, as MSstats'
#' `summaryMethod = "TMP"` does without imputation. Samples with no observed feature are
#' `NA`; a single-feature group returns that feature.
#'
#' Features observed in only some conditions bias additive summaries: median polish cannot
#' separate such a feature's offset from the condition effect, and the error, shared by
#' all samples of a condition, is invisible to the replicate variance. With
#' `conditions` given, features not observed in at least `min_per_condition` samples of
#' every condition are excluded before summarising (when a group has no feature left, all
#' its features are kept).
#'
#' @param M log2 matrix, features x samples.
#' @param group group (protein or site) of each row.
#' @param conditions optional condition of each column.
#' @param min_per_condition see above.
#' @return matrix, groups x samples, with attribute `n_features` (observed features used per
#'   group and sample).
#' @export
summarise_features <- function(M, group, conditions = NULL, min_per_condition = 1) {
  if (!is.null(conditions)) {
    conditions <- factor(conditions)
    ok <- Reduce(`&`, lapply(levels(conditions), function(l) rowSums(!is.na(M[, conditions == l, drop = FALSE])) >= min_per_condition))
    has_ok <- stats::ave(ok, group, FUN = any)
    M <- M[ok | !has_ok, , drop = FALSE]; group <- group[ok | !has_ok]
  }
  idx <- split(seq_len(nrow(M)), factor(group, levels = unique(group)))
  n_feat <- t(vapply(idx, function(i) colSums(!is.na(M[i, , drop = FALSE])), numeric(ncol(M))))
  out <- t(vapply(idx, function(i) {
    sub <- M[i, , drop = FALSE]
    sub <- sub[rowSums(!is.na(sub)) > 0, , drop = FALSE]
    if (!nrow(sub)) return(rep(NA_real_, ncol(M)))
    if (nrow(sub) == 1) return(sub[1, ])
    mp <- suppressWarnings(stats::medpolish(sub, na.rm = TRUE, trace.iter = FALSE, maxiter = 50))
    s <- mp$overall + mp$col
    s[colSums(!is.na(sub)) == 0] <- NA
    s
  }, numeric(ncol(M))))
  colnames(out) <- colnames(M); dimnames(n_feat) <- dimnames(out)
  attr(out, "n_features") <- n_feat
  out
}

#' Build a summarised `ptm_experiment` from simulated features
#'
#' @param sim output of [simulate_ptm_features()].
#' @param condition_filter exclude features not observed in every condition before
#'   summarising (see [summarise_features()]).
#' @return a `ptm_experiment` with the `truth` and `protein_truth` attributes expected by
#'   [evaluate_results()].
#' @export
sim_to_experiment <- function(sim, condition_filter = TRUE) {
  cnd <- if (condition_filter) sim$annotation$condition
  P <- summarise_features(sim$protein_features, sim$protein_map$protein, cnd)
  site_of <- sim$ptm_map$site
  Q <- summarise_features(sim$ptm_features, site_of, cnd)
  site_prot <- sim$ptm_map$protein[match(rownames(Q), site_of)]
  prot <- new_assay(P, data.frame(feature_id = rownames(P), protein_ids = rownames(P), gene = NA_character_))
  ptm <- new_assay(Q, data.frame(feature_id = rownames(Q), protein_ids = site_prot,
                                 stripped_sequence = NA_character_, modified_sequence = NA_character_,
                                 residue = "S", position = sub("^.*_S", "", rownames(Q)),
                                 localization_prob = 1, ptm_type = "Phospho", stringsAsFactors = FALSE))
  prot$n_features <- attr(P, "n_features"); ptm$n_features <- attr(Q, "n_features")
  x <- new_ptm_experiment(prot, ptm, annotation = list(protein = sim$annotation, ptm = sim$annotation),
                          params = list(simulation = sim$params))
  x <- add_log(x, paste0("Features summarised to proteins and sites by Tukey median polish",
                         if (condition_filter) " after excluding features not observed in every condition" else "", "."))
  attr(x, "truth") <- sim$truth
  attr(x, "protein_truth") <- sim$protein_truth
  x
}

#' MSstats long-format input from simulated features (for the native MSstatsPTM arm)
#'
#' @param sim output of [simulate_ptm_features()].
#' @return list(PTM = , PROTEIN = ) in MSstats label-free format, one row per feature x run,
#'   missing values as `NA` (MSstats' default `censoredInt = "NA"`). PTM `ProteinName`
#'   is `<protein>_<site>` so that MSstatsPTM finds the global protein.
#' @export
sim_to_msstats_input <- function(sim) {
  long <- function(M, prot_name, pep) {
    ann <- sim$annotation
    data.frame(ProteinName = rep(prot_name, times = ncol(M)), PeptideSequence = rep(pep, times = ncol(M)),
               PrecursorCharge = 2L, FragmentIon = NA_character_, ProductCharge = NA_integer_,
               IsotopeLabelType = "L", Condition = rep(ann$condition, each = nrow(M)),
               BioReplicate = rep(ann$sample, each = nrow(M)), Run = rep(ann$sample, each = nrow(M)),
               Fraction = 1L, Intensity = 2^as.vector(M), stringsAsFactors = FALSE)
  }
  list(PTM = long(sim$ptm_features, sim$ptm_map$site, sim$ptm_map$feature),
       PROTEIN = long(sim$protein_features, sim$protein_map$protein, sim$protein_map$feature))
}
