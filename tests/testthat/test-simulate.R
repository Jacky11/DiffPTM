test_that("imp4p missingness: rate per column and MCAR share are exact, MNAR targets low values", {
  set.seed(1)
  M <- matrix(rnorm(2000 * 4, 25, 2), 2000)
  for (pm in c(1, 0.5, 0.2)) {
    X <- add_missing_imp4p(M, rate = 0.3, pi_mcar = pm, para = 3)
    expect_equal(unname(colSums(is.na(X))), rep(600, 4))
    expect_equal(sum(attr(X, "mcar")), 4 * floor(600 * pm))
  }
  X <- add_missing_imp4p(M, rate = 0.3, pi_mcar = 0, para = 3)
  expect_lt(mean(M[is.na(X)]), mean(M) - 1)               # MNAR removes low intensities
  X <- add_missing_imp4p(M, rate = 0.3, pi_mcar = 1)
  expect_lt(abs(mean(M[is.na(X)]) - mean(M)), 0.2)         # MCAR does not
})

test_that("feature simulation: classes follow the applied effects and rates are respected", {
  sim <- simulate_ptm_features(n_proteins = 300, n_sites = 250, n_rep = 3, rate = c(0.1, 0.3), pi_mcar = 0.5, seed = 4)
  tr <- sim$truth
  expect_true(all(tr$dpu[tr$class %in% c("null", "protein_driven")] == 0))
  expect_true(all(tr$dpa[tr$class == "compensated"] == 0 & tr$dpu[tr$class == "compensated"] == 1))
  expect_equal(mean(is.na(sim$ptm_features)), 0.3, tolerance = 0.01)
  expect_equal(mean(is.na(sim$protein_features)), 0.1, tolerance = 0.01)
  expect_true(all(sim$ptm_map$protein %in% sim$protein_map$protein))
})

test_that("median polish summaries: one feature passes through, condition filter drops one-condition features", {
  M <- rbind(a = c(20, 21, 22, 23), b = c(25, 26, NA, NA), c = c(19, 20, 21, 22))
  s <- summarise_features(M[1, , drop = FALSE], "P1")
  expect_equal(unname(s[1, ]), M[1, ])
  s2 <- summarise_features(M, rep("P1", 3), conditions = c("A", "A", "B", "B"))
  expect_equal(unname(attr(s2, "n_features")[1, ]), c(2, 2, 2, 2))   # feature b excluded
  s3 <- summarise_features(M, rep("P1", 3))
  expect_equal(unname(attr(s3, "n_features")[1, ]), c(3, 3, 2, 2))
})

test_that("simulated experiment runs through the pipeline and the metrics are consistent", {
  sim <- simulate_ptm_features(n_proteins = 400, n_sites = 300, n_rep = 4, rho = 0.5, rate = c(0.05, 0.1), seed = 5)
  x <- match_ptm_to_protein(filter_min_observed(sim_to_experiment(sim), min_n = 2), "id")
  expect_false(is.null(x$ptm$protein_n_features))
  r <- method_limma(x, make_contrasts(x), pi0_method = "bh")
  ev <- evaluate_results(r, x, fdr = 0.05)
  dpu <- ev[ev$analysis == "ptm_DPU", ]
  expect_equal(dpu$TP + dpu$FP + dpu$FN + dpu$TN, dpu$n_features)
  expect_equal(dpu$FNR, 1 - dpu$sensitivity)
  expect_equal(dpu$accuracy, (dpu$TP + dpu$TN) / dpu$n_features)
  expect_true(dpu$auc > 0.5 && dpu$pauc <= 1 && dpu$pr_auc <= 1)
  expect_lt(dpu$observed_fdr, 0.15)
  bc <- evaluate_by_class(r, truth_table(x), 0.05)
  expect_true(all(c("null", "protein_driven", "usage_driven", "compensated") %in% bc$class[bc$analysis == "ptm_DPU"]))
  fc <- fdp_curve(r, truth_table(x))
  expect_true(all(diff(fc$n_called[fc$analysis == "ptm_DPU"]) >= 0))
  w <- method_limma(x, make_contrasts(x), analyses = "ptm_DPU", pi0_method = "bh", weights = TRUE)
  expect_equal(nrow(w), nrow(x$ptm$intensity))
})

test_that("MSstats inputs: native input names sites as protein_site and keeps NAs out", {
  sim <- simulate_ptm_features(n_proteins = 50, n_sites = 40, n_rep = 2, seed = 6)
  inp <- sim_to_msstats_input(sim)
  expect_true(all(sub("_S\\d+$", "", unique(inp$PTM$ProteinName)) %in% unique(inp$PROTEIN$ProteinName)))
  expect_true(all(c("ProteinName", "PeptideSequence", "Condition", "BioReplicate", "Run", "Intensity") %in% names(inp$PTM)))
  ann <- make_annotation(c("a1", "a2", "b1", "b2"), c("A", "A", "B", "B"), mixture = c("1", "2", "1", "2"))
  M <- matrix(c(1, NA, 3, 4), 1, 4, dimnames = list("x", ann$sample))
  lt <- DiffPTM:::msstats_long_tmt(M, "P1_S1", ann)
  expect_equal(nrow(lt), 3); expect_setequal(lt$Mixture, c("1", "2"))
  cm <- DiffPTM:::msstats_contrast_matrix(c("A", "B"), data.frame(name = "A_vs_B", group1 = "A", group2 = "B"))
  expect_equal(unname(cm[1, ]), c(1, -1))
})

test_that("mixture enters the limma design only without subject blocking", {
  d1 <- make_design(factor(c("A", "A", "B", "B")), mixture = c("1", "2", "1", "2"))
  expect_true(any(grepl("mixture", colnames(d1))))
  d2 <- make_design(factor(c("A", "A", "B", "B")), subject = c("s1", "s2", "s1", "s2"), mixture = c("1", "2", "1", "2"))
  expect_false(any(grepl("mixture", colnames(d2))))
})

test_that("simulate_ptm_features accepts per-site rho and a protein coefficient beta", {
  s1 <- simulate_ptm_features(n_rep = 3, rho = 0.5, seed = 1, n_proteins = 100, n_sites = 60)
  s2 <- simulate_ptm_features(n_rep = 3, rho = 0.5, seed = 1, n_proteins = 100, n_sites = 60, beta = 1)
  expect_identical(s1$ptm_features, s2$ptm_features)
  s3 <- simulate_ptm_features(n_rep = 3, rho = runif(60, 0, 0.9), seed = 1, n_proteins = 100, n_sites = 60)
  expect_equal(dim(s3$ptm_features), dim(s1$ptm_features))
  expect_error(simulate_ptm_features(n_rep = 3, rho = c(0.1, 0.2), n_proteins = 100, n_sites = 60), "rho")
  s4 <- simulate_ptm_features(n_rep = 3, rho = 0.5, seed = 1, n_proteins = 100, n_sites = 60, beta = 0.5)
  expect_false("compensated" %in% s4$truth$class)   # a compensated site no longer cancels its protein exactly
  expect_equal(s4$params$beta, 0.5)
})
