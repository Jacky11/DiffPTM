test_that("data model validates and reorders samples by annotation", {
  x <- simulate_experiment()
  expect_s3_class(x, "ptm_experiment")
  expect_equal(colnames(x$ptm$intensity), x$annotation$ptm$sample)
  expect_equal(levels(groups_of(x, "ptm")), c("A", "B"))
})

test_that("annotation guessing splits sample names", {
  ann <- guess_annotation(c("WT_t0_r1", "KO_t0_r1"), condition_pos = 1, factor_pos = 2, subject_pos = 3)
  expect_equal(ann$condition, c("WT", "KO")); expect_equal(ann$factor, c("t0", "t0")); expect_equal(ann$subject, c("r1", "r1"))
})

test_that("min-observed filter keeps features with enough values in one group", {
  x <- simulate_experiment(missing = 0.3)
  before <- nrow(x$ptm$intensity)
  y <- filter_min_observed(x, min_n = 3, assay = "ptm")
  obs <- !is.na(x$ptm$intensity)
  expected <- rowSums(obs[, 1:3]) >= 3 | rowSums(obs[, 4:6]) >= 3
  expect_equal(nrow(y$ptm$intensity), sum(expected))
  expect_lte(nrow(y$ptm$intensity), before)
  expect_true(any(grepl("kept", y$log)))
})

test_that("localization filter removes low-probability sites and keeps NA when asked", {
  x <- simulate_experiment()
  x$ptm$meta$localization_prob[1:3] <- 0.2; x$ptm$meta$localization_prob[4] <- NA
  y <- filter_localization(x, 0.75, keep_na = TRUE)
  expect_equal(nrow(y$ptm$intensity), 27)
  expect_equal(nrow(filter_localization(x, 0.75, keep_na = FALSE)$ptm$intensity), 26)
})

test_that("per-condition median normalisation aligns sample medians within groups", {
  x <- simulate_experiment()
  x$protein$intensity[, 2] <- x$protein$intensity[, 2] + 1
  y <- normalise(x, "protein", "median")
  med <- apply(y$protein$intensity[, 1:3], 2, median)
  expect_lt(diff(range(med)), 1e-8)
})

test_that("ID matching finds parents, averages shared, and drops unmatched", {
  x <- simulate_experiment()
  x$ptm$meta$protein_ids[1] <- "P1;P2"      # shared
  x$ptm$meta$protein_ids[2] <- "ZZZ"        # unmatched
  y <- match_ptm_to_protein(x, "id", shared = "mean")
  expect_equal(nrow(y$ptm$intensity), 30)                       # unmatched kept by default
  expect_true(all(is.na(y$ptm$protein_intensity[2, ])))
  expect_equal(nrow(match_ptm_to_protein(x, "id", drop_unmatched = TRUE)$ptm$intensity), 29)
  expect_equal(y$ptm$meta$n_matched_protein[1], 2L)
  expect_equal(unname(y$ptm$protein_intensity[1, 1]), mean(x$protein$intensity[1:2, 1]))
  expect_equal(nrow(match_ptm_to_protein(x, "id", shared = "drop", drop_unmatched = TRUE)$ptm$intensity), 28)
})

test_that("isoform suffixes are stripped when requested", {
  x <- simulate_experiment()
  x$ptm$meta$protein_ids[3] <- "P3-2"
  expect_equal(match_ptm_to_protein(x, "id")$ptm$meta$n_matched_protein[3], 1L)
  expect_equal(match_ptm_to_protein(x, "id", strip_isoforms = FALSE)$ptm$meta$n_matched_protein[3], 0L)
})

test_that("sequence matching finds peptides in a FASTA (I/L equivalence)", {
  fa <- tempfile(fileext = ".fasta")
  writeLines(c(">sp|P1|X", "MKKPEPTIDEAAA", ">sp|P2|Y", "MMMSEQLENCEKKK"), fa)
  seqs <- DiffPTM:::read_fasta(fa)
  hits <- DiffPTM:::find_peptides(c("PEPTIDE", "SEQIENCE", "NOTHERE"), seqs)
  expect_equal(hits[[1]], "P1"); expect_equal(hits[[2]], "P2"); expect_length(hits[[3]], 0)
})

test_that("collapse_to_sites sums peptidoforms of the same site", {
  x <- simulate_experiment()
  x$ptm$meta$protein_ids[2] <- "P1"; x$ptm$meta$position[2] <- "10"   # duplicate site of feature 1
  y <- collapse_to_sites(x, "sum")
  expect_equal(nrow(y$ptm$intensity), 29)
  expect_equal(y$ptm$meta$n_peptidoforms[1], 2L)
  expect_equal(unname(y$ptm$intensity[1, 1]), log2(2^x$ptm$intensity[1, 1] + 2^x$ptm$intensity[2, 1]))
})

test_that("detection discordance flags PTM observed without its protein", {
  x <- simulate_experiment()
  x$protein$intensity[1, 4:6] <- NA; x$protein$raw_detected[1, 4:6] <- FALSE
  y <- match_ptm_to_protein(x, "id")
  d <- detection_discordance(y, min_discordant = 2)
  expect_equal(d$flag[d$feature_id == "P1_S10" & d$group == "B"], "ptm_without_protein")
  expect_true(all(d$flag[d$feature_id != "P1_S10"] == "concordant"))
})

test_that("AhoCorasickTrie output (list of Keyword/Offset matches per text) is parsed", {
  hits <- list(list(list(Keyword = "PEPTLDE", Offset = 4L)), list(), list(list(Keyword = "PEPTLDE", Offset = 1L), list(Keyword = "SEQ", Offset = 2L)))
  r <- DiffPTM:::parse_ac_hits(hits, c("PEPTLDE", "SEQ", "NONE"), c("P1", "P2", "P3"))
  expect_equal(r$PEPTLDE, c("P1", "P3")); expect_equal(r$SEQ, "P3"); expect_null(r$NONE)
  df_hits <- list(data.frame(Keyword = "SEQ", Offset = 1L), data.frame(Keyword = character(0), Offset = integer(0)))
  expect_equal(DiffPTM:::parse_ac_hits(df_hits, "SEQ", c("A", "B"))$SEQ, "A")
  expect_null(DiffPTM:::parse_ac_hits(hits[1:2], "SEQ", c("A", "B", "C")))   # cannot align -> fallback
})

test_that("results carry degrees of freedom and a confidence interval consistent with the p-value", {
  sim <- simulate_ptm_features(n_rep = 3, rho = 0.5, rate = c(0.1, 0.2), seed = 2, n_proteins = 150, n_sites = 100)
  x <- match_ptm_to_protein(filter_min_observed(sim_to_experiment(sim), min_n = 2), "id")
  r <- run_methods(x, make_contrasts(x), methods = "limma", pi0_method = "bh")
  expect_true(all(c("df", "ci_lo", "ci_hi") %in% names(r)))
  u <- r[r$analysis == "ptm_DPU" & !is.na(r$pvalue), ]
  expect_true(all(u$df > 0))
  expect_equal((u$pvalue < 0.05), (u$ci_lo > 0 | u$ci_hi < 0))
})
