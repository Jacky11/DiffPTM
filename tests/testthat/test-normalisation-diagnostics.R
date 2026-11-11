mk <- function(seed = 1) {
  sim <- simulate_ptm_features(n_rep = 4, n_proteins = 600, n_sites = 800, seed = seed)
  filter_min_observed(sim_to_experiment(sim), min_n = 2)
}
test_that("no shift: global centering, case 1", {
  d <- normalisation_diagnostics(mk())
  expect_equal(d$recommendation, "global_median"); expect_equal(d$case, 1L)
  expect_setequal(d$table$assay, c("ptm", "protein"))
})
test_that("a technical shift of every site: global centering, case 3, proteome check reported", {
  x <- mk(); A <- groups_of(x, "ptm") == levels(groups_of(x, "ptm"))[1]
  x$ptm$intensity[, A] <- x$ptm$intensity[, A] + 0.6
  d <- normalisation_diagnostics(x)
  expect_equal(d$case, 3L); expect_equal(d$recommendation, "global_median")
  expect_match(d$message, "not shifted in the same way")
  x$protein$intensity[, A] <- x$protein$intensity[, A] + 0.6
  expect_match(normalisation_diagnostics(x)$message, "points to a loading difference")
})
test_that("a shift of a subset of sites, the bulk unchanged: per-condition centering, case 2", {
  x <- mk(2); A <- groups_of(x, "ptm") == levels(groups_of(x, "ptm"))[1]
  set.seed(3); sel <- sample(nrow(x$ptm$intensity), round(0.3 * nrow(x$ptm$intensity)))
  x$ptm$intensity[sel, A] <- x$ptm$intensity[sel, A] + 1.5
  d <- normalisation_diagnostics(x)
  expect_equal(d$case, 2L); expect_equal(d$recommendation, "median")
})
test_that("the diagnostic warns on normalised data", {
  expect_warning(normalisation_diagnostics(normalise(mk(), method = "global_median")), "already normalised")
})
