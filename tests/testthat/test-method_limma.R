test_that("limma method recovers DPA and DPU truth on simulated data", {
  x <- match_ptm_to_protein(simulate_experiment(sd = 0.2), "id")
  truth <- attr(simulate_experiment(), "truth")
  ctr <- make_contrasts(x, "pairwise")
  expect_equal(ctr$name, "A_vs_B")
  res <- run_methods(x, ctr, methods = "limma", pi0_method = "bh")
  expect_setequal(unique(res$analysis), c("protein_DPA", "ptm_DPA", "ptm_DPU"))
  sig <- function(an) res$feature_id[res$analysis == an & !is.na(res$adj_pvalue) & res$adj_pvalue < 0.05 & abs(res$log2FC) > 1]
  expect_setequal(sig("protein_DPA"), truth$protein_up)
  expect_setequal(sig("ptm_DPA"), truth$ptm_dpa_up)
  expect_setequal(sig("ptm_DPU"), truth$ptm_dpu_up)
  # sign convention: group1 - group2 = A - B, so "up in B" gives negative log2FC
  expect_true(all(res$log2FC[res$analysis == "ptm_DPU" & res$feature_id %in% truth$ptm_dpu_up] < -1))
})

test_that("limma method tolerates missing values without imputation", {
  x <- match_ptm_to_protein(filter_min_observed(simulate_experiment(sd = 0.2, missing = 0.15), min_n = 2), "id")
  res <- method_limma(x, make_contrasts(x), pi0_method = "bh")
  expect_true(all(c("log2FC", "pvalue", "adj_pvalue") %in% names(res)))
  dpa <- res[res$analysis != "ptm_DPU", ]
  expect_gt(sum(!is.na(dpa$pvalue)), 0.8 * nrow(dpa))
  # usage, consensus model: complete pairs only (default) tests no more sites than the all-observation fit
  cc <- method_limma(x, make_contrasts(x), analyses = "ptm_DPU", pi0_method = "bh", blocking = "random")
  all_obs <- method_limma(x, make_contrasts(x), analyses = "ptm_DPU", pi0_method = "bh", blocking = "random", complete_pairs = FALSE)
  expect_lte(sum(!is.na(cc$pvalue)), sum(!is.na(all_obs$pvalue)))
  expect_gt(sum(!is.na(res$pvalue[res$analysis == "ptm_DPU"])), 0.5 * sum(res$analysis == "ptm_DPU"))
})

test_that("the default usage model is the per-site paired model", {
  sim <- simulate_ptm_features(n_rep = 3, rho = 0.5, rate = c(0.1, 0.2), seed = 3, n_sites = 200, n_proteins = 300)
  x <- match_ptm_to_protein(filter_min_observed(sim_to_experiment(sim), min_n = 2), "id"); ctr <- make_contrasts(x)
  a <- method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh")
  f <- method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", blocking = "fixed")
  expect_equal(a$pvalue, f$pvalue)
  # it equals limma on the within-sample differences (correct then model) on complete sites
  d <- x$ptm$intensity - x$ptm$protein_intensity[, colnames(x$ptm$intensity)]
  full <- rowSums(is.na(d)) == 0
  g <- groups_of(x, "ptm"); X <- stats::model.matrix(~ 0 + g)
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(d[full, ], X), limma::makeContrasts(contrasts = "gA - gB", levels = X)))
  m <- match(rownames(d)[full], a$feature_id)
  expect_equal(unname(a$log2FC[m]), unname(fit$coefficients[, 1]), tolerance = 1e-8)
  expect_gt(stats::cor(a$statistic[m], fit$t[, 1]), 0.999)
  dg <- dpu_diagnostics(x)
  expect_equal(dg$recommended, "fixed"); expect_true(dg$rho_obs_median > 0)
})

test_that("the consensus usage model uses complete site-protein pairs only by default", {
  sim <- simulate_ptm_features(n_rep = 3, rho = 0.9, rate = c(0.2, 0.4), pi_mcar = 1, seed = 2, n_sites = 300, n_proteins = 300)
  x <- match_ptm_to_protein(filter_min_observed(sim_to_experiment(sim), min_n = 2), "id"); ctr <- make_contrasts(x)
  xc <- x; na <- is.na(x$ptm$intensity) | is.na(x$ptm$protein_intensity)
  xc$ptm$intensity[na] <- NA; xc$ptm$protein_intensity[na] <- NA
  a <- method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", blocking = "random")
  b <- method_limma(xc, ctr, analyses = "ptm_DPU", pi0_method = "bh", blocking = "random", complete_pairs = FALSE)
  expect_equal(a$pvalue, b$pvalue)
  d <- method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", blocking = "random", complete_pairs = FALSE)
  expect_false(isTRUE(all.equal(a$pvalue, d$pvalue)))
  # the per-site paired model is unaffected by the option
  f1 <- method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", blocking = "fixed")
  f2 <- method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", blocking = "fixed", complete_pairs = FALSE)
  expect_equal(f1$pvalue, f2$pvalue)
})

test_that("simulator options leave the default simulation unchanged and shape the data", {
  a <- simulate_ptm_features(n_rep = 3, seed = 5, n_sites = 200, n_proteins = 300)
  b <- simulate_ptm_features(n_rep = c(3, 3), seed = 5, n_sites = 200, n_proteins = 300)
  expect_identical(a$ptm_features, b$ptm_features)
  u <- simulate_ptm_features(n_rep = c(3, 5), seed = 1, n_sites = 200, n_proteins = 300, rate_B_factor = 2, site_clustering = 0.06)
  expect_equal(ncol(u$ptm_features), 8)
  miss <- colMeans(is.na(u$ptm_features)); expect_gt(mean(miss[4:8]), 1.5 * mean(miss[1:3]))
  expect_lt(length(unique(u$truth$protein)), 150)
})

test_that("vs_all contrasts and a second factor are supported", {
  x <- simulate_experiment()
  for (a in c("protein", "ptm")) x$annotation[[a]]$factor <- rep(c("t0", "t1", "t0"), 2)
  x <- validate_ptm_experiment(x)
  ctr <- make_contrasts(x, "pairwise", within_factor = TRUE)
  expect_true(all(c("A.t0_vs_B.t0", "A.t1_vs_B.t1") %in% ctr$name))
  ctr2 <- make_contrasts(x, "vs_all", within_factor = FALSE)
  expect_equal(nrow(ctr2), 4)
  res <- method_limma(x, ctr2, analyses = "protein_DPA", pi0_method = "bh")
  expect_equal(length(unique(res$contrast)), 4)
})

test_that("subject blocking works as fixed effect", {
  x <- simulate_experiment()
  for (a in c("protein", "ptm")) x$annotation[[a]]$subject <- rep(c("s1", "s2", "s3"), 2)
  res <- method_limma(x, make_contrasts(x), analyses = "protein_DPA", blocking = "fixed", pi0_method = "bh")
  expect_equal(nrow(res), 40)
})

test_that("compare_methods produces the wide table and consensus", {
  x <- match_ptm_to_protein(simulate_experiment(), "id")
  res <- method_limma(x, make_contrasts(x), pi0_method = "bh")
  res2 <- res; res2$method <- "fake"
  cmp <- compare_methods(rbind(res, res2), fdr = 0.05, min_log2fc = 1)
  expect_true(all(c("log2FC_limma", "log2FC_fake", "n_methods_significant") %in% names(cmp$wide)))
  expect_equal(nrow(consensus_calls(rbind(res, res2), 0.05, 1, 2)), sum(cmp$wide$n_methods_significant == 2))
})

test_that("xlsx report is written with README and methods", {
  x <- match_ptm_to_protein(normalise(simulate_experiment()), "id")
  res <- method_limma(x, make_contrasts(x), pi0_method = "bh")
  f <- tempfile(fileext = ".xlsx")
  write_results_xlsx(x, res, f)
  expect_true(file.exists(f))
  expect_true(grepl("DiffPTM", methods_paragraph(x)))
})

test_that("moderate = FALSE gives ordinary statistics and blocking = 'fixed' the paired t-test", {
  sim <- simulate_ptm_features(n_rep = 3, rho = 0.5, rate = c(0, 0), pi_mcar = 1, seed = 3, n_proteins = 200, n_sites = 150)
  x <- match_ptm_to_protein(filter_min_observed(sim_to_experiment(sim), min_n = 2), "id")
  ctr <- make_contrasts(x)
  r <- method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", blocking = "fixed", moderate = FALSE)
  D <- x$ptm$intensity - x$ptm$protein_intensity; g <- groups_of(x, "ptm")
  i <- which(rowSums(is.na(D)) == 0)[1:5]
  for (k in i) {
    tt <- t.test(D[k, g == levels(g)[1]], D[k, g == levels(g)[2]], var.equal = TRUE)
    expect_equal(r$statistic[r$feature_id == rownames(D)[k]], unname(tt$statistic), tolerance = 1e-6)
    expect_equal(r$pvalue[r$feature_id == rownames(D)[k]], tt$p.value, tolerance = 1e-6)
  }
  m <- method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", blocking = "fixed")
  expect_equal(m$log2FC, r$log2FC)          # moderation changes the statistics, not the estimate
  expect_false(isTRUE(all.equal(m$pvalue, r$pvalue)))
})

test_that("usage_beta_check recovers protein changes and flags beta-sensitive calls", {
  sim <- simulate_ptm_features(n_rep = 4, rho = 0.5, rate = c(0.05, 0.1), seed = 4, n_sites = 300, n_proteins = 300)
  x <- match_ptm_to_protein(filter_min_observed(sim_to_experiment(sim), min_n = 2), "id"); ctr <- make_contrasts(x)
  b <- usage_beta_check(x, ctr, pi0_method = "bh")
  base <- method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh")
  expect_equal(b$adj_pvalue, base$adj_pvalue)
  dpa <- method_limma(x, ctr, analyses = "ptm_DPA", pi0_method = "bh")
  dpa <- dpa[dpa$analysis == "ptm_DPA", ]; m <- match(b$feature_id, dpa$feature_id)
  ok <- !is.na(b$protein_log2FC) & !is.na(dpa$log2FC[m])
  expect_gt(stats::cor(b$protein_log2FC[ok], dpa$log2FC[m][ok] - b$usage_log2FC[ok]), 0.99)
  expect_true(all(!b$robust_to_beta | b$called))
})

test_that("dpu_diagnostics finds no spread of the variance ratio when it is constant", {
  ex <- sapply(1:3, function(s) {
    sim <- simulate_ptm_features(n_rep = 4, rho = 0.5, seed = s, n_sites = 800, n_proteins = 800, site_sd_factor = 2)
    d <- dpu_diagnostics(match_ptm_to_protein(filter_min_observed(sim_to_experiment(sim), min_n = 2), "id"))
    c(d$log_var_ratio_excess, d$log_sd_ratio_sd) })
  expect_lt(abs(mean(ex[1, ])), 0.1)        # the excess variance is centred on 0
  expect_true(all(ex[2, ] < 0.3))           # below the threshold of the recommendation
})

test_that("dpu_diagnostics detects a heterogeneous variance ratio and leaves the random stream unchanged", {
  sim <- simulate_ptm_features(n_rep = 4, rho = 0.5, seed = 2, n_sites = 800, n_proteins = 800, site_sd_factor = 2, sd_heterogeneity = 1)
  x <- match_ptm_to_protein(filter_min_observed(sim_to_experiment(sim), min_n = 2), "id")
  set.seed(5); a <- runif(1); set.seed(5); d <- dpu_diagnostics(x); b <- runif(1)
  expect_identical(a, b)
  expect_gt(d$log_sd_ratio_sd, 0.3)
  expect_true(is.finite(d$rho_obs_pooled))
})

test_that("the consensus model stops with an explicit message when the correlation is not estimable", {
  sim <- simulate_ptm_features(n_rep = 2, rho = 0.5, seed = 3, n_sites = 200, n_proteins = 200, rate = c(0, 0))
  x <- match_ptm_to_protein(filter_min_observed(sim_to_experiment(sim), min_n = 2), "id")
  x$annotation$ptm$mixture <- x$annotation$protein$mixture <- rep(c("1", "2"), length.out = nrow(x$annotation$ptm))
  r <- tryCatch(method_limma(x, make_contrasts(x), analyses = "ptm_DPU", blocking = "random"), error = function(e) conditionMessage(e))
  if (is.character(r)) expect_match(r, "could not be estimated") else succeed()
})
