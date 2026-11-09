#' Diagnostics for the choice of pairing in the usage (DPU) model
#'
#' The per-site paired single model (default of [method_limma()]) analyses each site's
#' within-sample site-minus-protein differences and makes no assumption on the variances.
#' The consensus single model (`blocking = "random"`) pools the residual variances of the
#' site and of its protein into one, and shares one within-sample correlation between all
#' sites; it doubles the residual degrees of freedom, which helps with two replicates per
#' group, but it is anti-conservative when the site and its protein have different, and
#' differently different, residual variances. This function reports the quantities that
#' decide between the two, computed on the complete site-protein pairs after centring each
#' value within its group (and TMT mixture):
#'
#' * `rho_obs_median`, `rho_obs_q25`, `rho_obs_q75`: per-site Pearson correlation of the
#'   site and protein residuals (sites with at least `min_pairs` complete pairs);
#' * `log_var_ratio_median`: median over sites of `log(s_S^2 / s_P^2)` (natural log; 0 when the
#'   site and its protein vary alike, about 1.6-2.4 in the phosphoproteomes of the DiffPTM
#'   benchmark);
#' * `log_sd_ratio_sd`: the spread over sites of the true log SD ratio,
#'   `sqrt(max(V_obs - V_0, 0)) / 2`, where `V_obs` is the variance over sites of the observed
#'   `log(s_S^2 / s_P^2)` and `V_0` its expectation when every site has the same variance ratio.
#'   For a site with `nu` residual degrees of freedom and residual correlation `r` between site
#'   and protein, the sampling variance of the log ratio is `2 * trigamma(nu / 2)` when `r = 0`
#'   and smaller when `r != 0` (the two variances then share sampling error);
#'   `V_0` is the mean over sites of this variance, computed by Monte Carlo (`n_mc` draws of
#'   two correlated chi-squared sums with `nu` degrees of freedom and correlation `rho_obs_pooled`,
#'   with a fixed internal seed that leaves the user's random stream unchanged). The factor 1/2
#'   converts the log variance ratio to the log SD ratio;
#' * `log_var_ratio_excess`: `V_obs - V_0` (close to 0, positive or negative, when the ratio is the
#'   same for every site; checked in the DiffPTM benchmark, `review5_diagnostic_null.R`);
#' * `rho_obs_pooled`: correlation of the site and protein residuals pooled over sites;
#' * `replicates_per_group`: the median number of complete pairs per condition;
#' * `prior_df_per_site`: the prior degrees of freedom `d0` of the moderated variances of the
#'   within-sample differences (`limma::squeezeVar`);
#' * `recommended`: `"fixed"` (per-site paired) unless there are fewer than three replicates
#'   per group and the site and protein variances are similar enough for pooling
#'   (`abs(log_var_ratio_median) <= max_log_var_ratio` and `log_sd_ratio_sd <= max_log_sd_ratio_sd`),
#'   in which case `"random"` (consensus) gains degrees of freedom at little risk.
#'
#' The thresholds come from the simulations of the DiffPTM benchmark: the consensus model was
#' calibrated when the variance ratio was homogeneous over sites (estimated spread 0.01-0.12)
#' and anti-conservative when it varied as in the real datasets (spread 0.4-0.75; real data 0.54-0.99).
#'
#' @param x a `ptm_experiment` after [match_ptm_to_protein()].
#' @param min_pairs minimal number of complete pairs for the per-site correlation.
#' @param max_log_var_ratio,max_log_sd_ratio_sd thresholds of the recommendation.
#' @param n_mc Monte Carlo draws per distinct number of residual degrees of freedom for `V_0`.
#' @return a one-row data frame.
#' @examples
#' sim <- simulate_ptm_features(n_rep = 3, n_proteins = 300, n_sites = 200, seed = 1)
#' x <- match_ptm_to_protein(filter_min_observed(sim_to_experiment(sim), min_n = 2), "id")
#' dpu_diagnostics(x)
#' @export
dpu_diagnostics <- function(x, min_pairs = 6, max_log_var_ratio = log(3), max_log_sd_ratio_sd = 0.3, n_mc = 2000) {
  if (is.null(x$ptm$protein_intensity)) stop("Run match_ptm_to_protein() first")
  S <- x$ptm$intensity; P <- x$ptm$protein_intensity[, colnames(S), drop = FALSE]
  a <- x$annotation$ptm
  grp <- paste(as.character(groups_of(x, "ptm")), if (!is.null(a$mixture) && !all(is.na(a$mixture))) a$mixture else "")
  ok <- !is.na(S) & !is.na(P); S[!ok] <- NA; P[!ok] <- NA
  cen <- function(M) { for (g in unique(grp)) { j <- grp == g; M[, j] <- M[, j] - rowMeans(M[, j, drop = FALSE], na.rm = TRUE) }; M }
  rS <- cen(S); rP <- cen(P)
  cond <- as.character(groups_of(x, "ptm"))
  count_by <- function(f) { m <- vapply(unique(f), function(g) rowSums(ok[, f == g, drop = FALSE]), numeric(nrow(S)))
    if (is.null(dim(m))) matrix(m, nrow = nrow(S)) else m }
  per_grp <- count_by(cond)                          # complete pairs per condition
  nu <- rowSums(ok) - rowSums(count_by(grp) > 0)     # residual df after condition (x mixture) means
  keep <- apply(per_grp, 1, min) >= 2 & nu >= 2
  vS <- rowSums(rS^2, na.rm = TRUE) / nu; vP <- rowSums(rP^2, na.rm = TRUE) / nu
  k2 <- keep & vS > 0 & vP > 0
  lr <- log(vS[k2] / vP[k2])
  rho_pool <- stats::cor(as.vector(rS[k2, ]), as.vector(rP[k2, ]), use = "complete.obs")
  v0 <- log_ratio_null_var(nu[k2], rho_pool, n_mc)
  excess <- stats::var(lr) - v0
  rr <- vapply(which(keep & rowSums(ok) >= min_pairs), function(i) stats::cor(rS[i, ], rP[i, ], use = "complete.obs"), numeric(1))
  vd <- rowSums((rS - rP)^2, na.rm = TRUE) / nu
  d0 <- tryCatch(limma::squeezeVar(vd[keep & vd > 0], nu[keep & vd > 0])$df.prior, error = function(e) NA_real_)
  reps <- stats::median(per_grp[keep, , drop = FALSE])
  sd_ratio <- sqrt(max(excess, 0)) / 2
  rec <- if (reps < 3 && abs(stats::median(lr)) <= max_log_var_ratio && sd_ratio <= max_log_sd_ratio_sd) "random" else "fixed"
  data.frame(sites = sum(k2), replicates_per_group = reps,
             rho_obs_median = stats::median(rr), rho_obs_q25 = unname(stats::quantile(rr, 0.25)), rho_obs_q75 = unname(stats::quantile(rr, 0.75)), rho_obs_pooled = rho_pool,
             log_var_ratio_median = stats::median(lr), log_var_ratio_excess = excess, log_sd_ratio_sd = sd_ratio, prior_df_per_site = d0,
             recommended = rec, stringsAsFactors = FALSE)
}

# Expected variance over sites of log(s_S^2 / s_P^2) when the variance ratio is the same for every site:
# mean over sites of the sampling variance for nu residual df and residual correlation rho (Monte Carlo,
# fixed internal seed; the global random stream is restored).
log_ratio_null_var <- function(nu, rho, n_mc = 2000) {
  if (!length(nu)) return(NA_real_)
  if (!is.finite(rho)) rho <- 0
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv()))
  set.seed(20261108)
  u <- sort(unique(nu)); w <- tabulate(match(nu, u), length(u)) / length(nu)
  v <- vapply(u, function(n) {
    z1 <- matrix(stats::rnorm(n_mc * n), n_mc); z2 <- rho * z1 + sqrt(1 - rho^2) * matrix(stats::rnorm(n_mc * n), n_mc)
    stats::var(log(rowSums(z1^2) / rowSums(z2^2)))
  }, numeric(1))
  sum(w * v)
}
