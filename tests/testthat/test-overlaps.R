x0 <- function() { sim <- simulate_ptm_features(n_rep = 4, n_proteins = 300, n_sites = 250, seed = 5); filter_min_observed(sim_to_experiment(sim), min_n = 2) }
test_that("Venn and UpSet plots are drawn from sets", {
  set.seed(1); S <- matrix(runif(400) > 0.3, 100, 4, dimnames = list(NULL, letters[1:4]))
  for (k in 2:4) expect_s3_class(plot_venn(S[, 1:k]), "ggplot")
  expect_s3_class(plot_upset(S), "ggplot")
  # the numbers drawn in a Venn diagram add up to the elements in at least one set
  b <- ggplot2::ggplot_build(plot_venn(S[, 1:3]))
  txt <- b$data[[which(vapply(b$data, function(d) "label" %in% names(d), TRUE))[1]]]
  expect_equal(sum(as.integer(txt$label)), sum(rowSums(S[, 1:3]) > 0))
})
test_that("detection overlaps, dendrograms, scree and condition means are drawn", {
  x <- x0()
  expect_s3_class(plot_detection_overlap(x, "ptm"), "ggplot")
  expect_s3_class(plot_detection_overlap(x, "ptm", group = "A"), "ggplot")            # 4 samples: Venn
  expect_s3_class(plot_detection_overlap(x, "ptm", level = "protein"), "ggplot")
  expect_lte(nrow(detection_sets(x, "ptm", "protein")), nrow(detection_sets(x, "ptm")))
  for (m in c("pearson", "jaccard")) expect_s3_class(plot_dendrogram(x, "protein", m), "ggplot")
  expect_s3_class(plot_pca_scree(x, "protein"), "ggplot")
  expect_match(plot_jaccard_matrix(x, "ptm")$labels$subtitle, "lowest")
  x <- match_ptm_to_protein(normalise(x), "id"); ctr <- make_contrasts(x)
  r <- suppressMessages(run_methods(x, ctr, methods = "limma"))
  expect_s3_class(plot_condition_means(x, r[r$analysis == "ptm_DPA", ], ctr$group1[1], ctr$group2[1], "ptm"), "ggplot")
  f <- tempfile(fileext = ".xlsx"); write_results_xlsx(x, r, f)
  d <- openxlsx::read.xlsx(f, "ptm_DPU_limma")
  expect_true(all(c("mean_log2_A", "mean_log2_B", "protein_mean_log2_A") %in% names(d)))
})
