test_that("simulation carries consistent truth and the limma method recovers it", {
  x <- simulate_ptm_experiment(n_proteins = 300, n_ptm = 200, n_rep = 4, missing = c(0, 0), seed = 2)
  tr <- attr(x, "truth")
  expect_equal(nrow(tr), 200)
  expect_setequal(unique(tr$class), c("null", "protein_driven", "usage_driven", "compensated"))
  expect_true(all(tr$dpu[tr$class == "protein_driven"] == 0))
  expect_true(all(tr$dpa[tr$class == "compensated"] == 0 & tr$dpu[tr$class == "compensated"] == 1))
  x <- match_ptm_to_protein(x, "id")
  res <- run_methods(x, make_contrasts(x), "limma", pi0_method = "bh")
  ev <- evaluate_results(res, x, fdr = 0.05)
  dpu <- ev[ev$analysis == "ptm_DPU", ]
  expect_lt(dpu$observed_fdr, 0.15)
  expect_gt(dpu$sensitivity, 0.6)
  expect_gt(dpu$auc, 0.95)
  expect_lt(abs(dpu$fc_bias), 0.2)
})

test_that("missing values are generated at the requested rate with MNAR bias", {
  x <- simulate_ptm_experiment(n_proteins = 300, n_ptm = 300, n_rep = 3, missing = c(0.1, 0.3), mnar = 1, seed = 3)
  expect_equal(mean(is.na(x$ptm$intensity)), 0.3, tolerance = 0.02)
  expect_equal(mean(is.na(x$protein$intensity)), 0.1, tolerance = 0.02)
  Q <- x$ptm$intensity; obs <- Q[!is.na(Q)]
  # MNAR: the observed values must be higher on average than the full matrix would be (base - 3 ~ 21)
  expect_gt(mean(obs), 21)
})
