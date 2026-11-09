sim_pair <- function(seed = 2) {
  sim <- simulate_ptm_features(n_rep = 3, n_proteins = 150, n_sites = 100, seed = seed)
  x <- sim_to_experiment(sim)
  x$protein$n_features <- NULL; x$ptm$n_features <- NULL
  x
}
add_empty_runs <- function(a, runs) {          # runs listed in a table but without any value
  I <- cbind(a$intensity, matrix(NA_real_, nrow(a$intensity), length(runs), dimnames = list(NULL, runs)))
  a$intensity <- I; a$raw_detected <- !is.na(I); a
}
used <- function(d) d[nzchar(d$condition), ]

test_that("combined MaxQuant-like search: empty runs dropped, prefix removed, runs paired", {
  x <- sim_pair(); s <- colnames(x$ptm$intensity)
  ptm <- add_empty_runs(set_run_names(x$ptm, paste0("Phospho_", s)), c(s, paste0("Ubi_", s)))
  g <- suggest_design(x$protein, ptm)
  u <- used(g$annotation$ptm)
  expect_setequal(u$run, paste0("Phospho_", s))
  expect_setequal(u$sample, s)
  expect_equal(length(unique(u$condition)), 2)
  expect_setequal(used(g$annotation$protein)$sample, s)
})

test_that("short fraction markers (P_ / T_) without any keyword are removed and runs paired", {
  x <- sim_pair(); s <- colnames(x$ptm$intensity)
  g <- suggest_design(set_run_names(x$protein, paste0("T_", s)), set_run_names(x$ptm, paste0("P_", s)))
  expect_setequal(used(g$annotation$ptm)$sample, s)
  expect_setequal(used(g$annotation$protein)$sample, s)
})

test_that("a combined search with non-empty total runs in the PTM table is split by the number of sites", {
  x <- sim_pair(); s <- colnames(x$ptm$intensity)
  enr <- paste0("E_", s); tot <- paste0("T_", s)
  ptm <- x$ptm; I <- ptm$intensity
  sparse <- I; sparse[sample(length(sparse), round(0.97 * length(sparse)))] <- NA     # total runs: few sites
  ptm$intensity <- cbind(`colnames<-`(I, enr), `colnames<-`(sparse, tot)); ptm$raw_detected <- !is.na(ptm$intensity)
  prot <- x$protein; P <- prot$intensity
  prot$intensity <- cbind(`colnames<-`(P, enr), `colnames<-`(P, tot)); prot$raw_detected <- !is.na(prot$intensity)
  g <- suggest_design(prot, ptm)
  expect_setequal(used(g$annotation$ptm)$run, enr)
  expect_setequal(used(g$annotation$protein)$run, tot)
  expect_setequal(intersect(used(g$annotation$ptm)$sample, used(g$annotation$protein)$sample), s)
  expect_true(any(grepl("Combined search", g$messages)))
})

test_that("DIA-NN-like file names: fraction tokens and dates differ, runs still paired", {
  x <- sim_pair(); s <- colnames(x$ptm$intensity)                    # A_1, A_2, ... B_3
  prot <- set_run_names(x$protein, paste0("20230315_QE_FP_", s))
  ptm <- set_run_names(x$ptm, paste0("20230317_QE_Phos_", s))
  g <- suggest_design(prot, ptm)
  up <- used(g$annotation$protein); uq <- used(g$annotation$ptm)
  expect_equal(length(intersect(up$sample, uq$sample)), length(s))
  expect_equal(uq$sample[match(paste0("20230317_QE_Phos_", s[1]), uq$run)], up$sample[match(paste0("20230315_QE_FP_", s[1]), up$run)])
})

test_that("a re-injection label is paired with the closest proteome run", {
  x <- sim_pair(3); s <- colnames(x$ptm$intensity)
  new <- paste0("Phospho-", s); new[2] <- paste0(new[2], "_Octobre")
  g <- suggest_design(set_run_names(x$protein, paste0("Proteome-", s)), set_run_names(x$ptm, new))
  uq <- used(g$annotation$ptm)
  expect_setequal(uq$sample, used(g$annotation$protein)$sample)
})

test_that("conditions are read wherever the replicate number is", {
  expect_equal(conditions_from_names(c("R1-A", "R1-B", "R2-A", "R2-B")), c("A", "B", "A", "B"))
  expect_equal(conditions_from_names(c("mock_4h_1", "mock_4h_2", "inf_24h_1", "inf_24h_2")), c("mock_4h", "mock_4h", "inf_24h", "inf_24h"))
  expect_equal(conditions_from_names(c("mock_24h_1", "mock_24h_2", "inf_24h_1", "inf_24h_2")), c("mock", "mock", "inf", "inf"))  # parts common to all names dropped
  expect_equal(conditions_from_names(c("01_NS_R1", "02_TNF_R1", "03_NS_R2", "04_TNF_R2")), c("NS", "TNF", "NS", "TNF"))
  expect_equal(conditions_from_names(c("L_3", "H_3", "L_4", "H_4")), c("L", "H", "L", "H"))
  expect_true(all(is.na(conditions_from_names(c("a", "b", "c")))))
})

test_that("MaxQuant SILAC columns become label x experiment samples", {
  m <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), 1, 8, dimnames = list(NULL, c("L", "H", "3", "L 3", "H 3", "4", "L 4", "H 4")))
  expect_equal(colnames(maxquant_label_samples(m)), c("L_3", "H_3", "L_4", "H_4"))
  m[, c("L 3", "H 3", "L 4", "H 4")] <- 0
  expect_equal(colnames(maxquant_label_samples(m)), c("L", "H"))
  lf <- matrix(1, 1, 2, dimnames = list(NULL, c("ctrl_1", "ctrl_2")))
  expect_identical(maxquant_label_samples(lf), lf)
})

test_that("suggest_design leaves matching names alone (example dataset)", {
  ex <- system.file("extdata", "example", package = "DiffPTM")
  g <- suggest_design(import_protein(file.path(ex, "protein.tsv"), "Generic"), import_ptm(file.path(ex, "ptm.tsv"), "Generic"))
  expect_identical(g$annotation$ptm$run, g$annotation$ptm$sample)
  expect_setequal(unique(g$annotation$ptm$condition), c("control", "treated"))
})

test_that("the design module proposes a design and applies it", {
  x <- sim_pair(); s <- colnames(x$ptm$intensity)
  ptm <- add_empty_runs(set_run_names(x$ptm, paste0("Phospho_", s)), s)
  raw <- shiny::reactiveVal(list(protein = x$protein, ptm = ptm))
  shiny::testServer(mod_design_server, args = list(raw = raw), {
    session$flushReact()
    expect_setequal(used(ann$ptm)$sample, s)
    session$setInputs(apply = 1)
    e <- session$getReturned()()
    expect_s3_class(e, "ptm_experiment")
    expect_setequal(colnames(e$ptm$intensity), s)
    expect_setequal(intersect(e$annotation$ptm$sample, e$annotation$protein$sample), s)
  })
})
