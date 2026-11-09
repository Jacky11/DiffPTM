test_that("merge_redundant_ptm merges identical rows reported for several proteins", {
  I <- rbind(c(10, 11, NA, 12), c(10, 11, NA, 12), c(10, 11, NA, 12), c(9, 9, 9, 9), c(10, 11, 13, 12))
  meta <- data.frame(feature_id = c("A0A087WU74_T648", "Q96NU1_T485", "H7BY14_T393", "P12345_S10", "P99999_S5"),
                     protein_ids = c("A0A087WU74", "Q96NU1", "H7BY14", "P12345", "P99999"),
                     modified_sequence = c("AVGCR", "AVGCR", "AVGCR", "KLMSR", "AVGCR"), multiplicity = 1)
  colnames(I) <- paste0("s", 1:4)
  ann <- make_annotation(colnames(I), c("a", "a", "b", "b"))
  x <- new_ptm_experiment(ptm = new_assay(I, meta), annotation = list(ptm = ann))
  y <- merge_redundant_ptm(x)
  expect_equal(nrow(y$ptm$intensity), 3)
  expect_true("Q96NU1_T485" %in% y$ptm$meta$feature_id)          # reviewed accession kept
  expect_equal(y$ptm$meta$n_redundant[y$ptm$meta$feature_id == "Q96NU1_T485"], 3)
  expect_setequal(strsplit(y$ptm$meta$protein_ids[y$ptm$meta$feature_id == "Q96NU1_T485"], ";")[[1]],
                  c("A0A087WU74", "Q96NU1", "H7BY14"))
  expect_true("P99999_S5" %in% y$ptm$meta$feature_id)            # same sequence, other values: kept
})
