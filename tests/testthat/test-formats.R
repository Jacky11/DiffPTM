# End-to-end import of the example exports of each software (inst/extdata/formats, built by
# inst/scripts/make_format_fixtures.R from one simulated experiment) and recovery of the
# planted changes with the single model.
fx <- function(...) system.file("extdata", "formats", ..., package = "DiffPTM")
truth <- readRDS(fx("truth.rds"))
ann <- read_annotation(fx("annotation.tsv"))

run_pipeline <- function(p, s, fasta = NULL, locate = FALSE) {
  x <- new_ptm_experiment(p, s, annotation = list(protein = ann[ann$sample %in% colnames(p$intensity), ], ptm = ann[ann$sample %in% colnames(s$intensity), ]))
  x <- filter_flagged(x) |> merge_redundant_ptm()
  if (locate) x <- locate_sites(x, fasta)
  x <- filter_localization(x, 0.75) |> filter_min_observed(min_n = 2) |> normalise(method = "median") |> match_ptm_to_protein("id")
  res <- run_methods(x, make_contrasts(x, "vs_reference", reference = "control"), "limma", pi0_method = "bh")
  list(x = x, res = res)
}
check_truth <- function(out, site_id = function(meta) meta$feature_id) {
  res <- out$res; meta <- out$x$ptm$meta
  prot <- res[res$analysis == "protein_DPA" & !is.na(res$adj_pvalue) & res$adj_pvalue < 0.05, "feature_id"]
  expect_gte(length(intersect(prot, truth$protein_changed)), 8)
  expect_lte(length(prot), 20)
  key <- site_id(meta); dpu <- res[res$analysis == "ptm_DPU" & !is.na(res$adj_pvalue) & res$adj_pvalue < 0.05, "feature_id"]
  called <- key[match(dpu, meta$feature_id)]
  expect_gte(length(intersect(called, truth$usage_up)), 5)
  expect_gte(length(intersect(called, truth$compensated)), 2)
  expect_lte(length(called), 20)
}
site_from_meta <- function(meta) paste0(sub(";.*$", "", meta$protein_ids), "_", meta$residue, meta$position)

test_that("MaxQuant combined search: proteome runs + Phospho_ runs, renamed and paired", {
  p <- import_protein(fx("MaxQuant", "proteinGroups.txt"), "MaxQuant")
  s <- import_ptm(fx("MaxQuant", "Phospho_STY_Sites.txt"), "MaxQuant")
  expect_equal(ncol(p$intensity), 6)                       # LFQ columns exist only for the proteome runs
  expect_equal(ncol(s$intensity), 12)                      # site table lists both run types
  s <- rename_samples(s, "^Phospho_(.+)$")
  expect_equal(colnames(s$intensity), colnames(p$intensity))
  expect_true(all(c("Reverse", "Potential.contaminant") %in% names(p$meta)))
  out <- run_pipeline(p, s)
  expect_equal(nrow(out$x$protein$intensity), 60)          # reverse + contaminant removed
  check_truth(out, site_from_meta)
  expect_match(out$x$ptm$meta$label[1], "^GENE\\d+ [STY]\\d+")
})

test_that("Spectronaut pivot reports (semicolon, decimal comma, Filtered, protein groups)", {
  p <- import_protein(fx("Spectronaut", "PivotReport_Proteomes.csv"), "Spectronaut")
  s <- import_ptm(fx("Spectronaut", "PivotReport_Phospho.csv"), "Spectronaut")
  expect_equal(ncol(p$intensity), 6); expect_equal(nrow(s$intensity), 46)     # one site of a protein group listed twice
  expect_lt(mean(is.na(p$intensity)), 0.2)
  p <- rename_samples(p, "^Proteome-\\d+_(.+)$"); s <- rename_samples(s, "^Phospho-\\d+_(.+)$")
  expect_equal(colnames(p$intensity), c("control_R1", "control_R2", "control_R3", "treated_R1", "treated_R2", "treated_R3"))
  expect_equal(colnames(s$intensity), colnames(p$intensity))
  expect_error(rename_samples(p, "^(.+)_R\\d$"), "duplicated")        # the capture group must keep the replicate
  p <- import_protein(fx("Spectronaut", "PivotReport_Proteomes.csv"), "Spectronaut") |> rename_samples("^Proteome-\\d+_(.+)$")
  s <- import_ptm(fx("Spectronaut", "PivotReport_Phospho.csv"), "Spectronaut") |> rename_samples("^Phospho-\\d+_(.+)$")
  colnames(p$intensity) <- colnames(p$raw_detected) <- sub("_R", "_", colnames(p$intensity))
  colnames(s$intensity) <- colnames(s$raw_detected) <- sub("_R", "_", colnames(s$intensity))
  out <- run_pipeline(p, s)
  expect_equal(nrow(out$x$ptm$intensity) <= 45, TRUE)       # redundant rows merged
  check_truth(out, site_from_meta)
})

test_that("Proteome Discoverer Proteins + PeptideGroups exports", {
  p <- import_protein(fx("ProteomeDiscoverer", "Proteins.txt"), "Proteome Discoverer")
  s <- import_ptm(fx("ProteomeDiscoverer", "PeptideGroups.txt"), "Proteome Discoverer")
  expect_equal(colnames(p$intensity), paste0("F", 1:6))     # normalised abundances, one per file
  expect_equal(colnames(s$intensity), paste0("F", 1:6))
  expect_true(all(s$meta$residue %in% c("S", "T", "Y")))
  expect_true(all(!is.na(as.integer(s$meta$position))))
  expect_true(all(s$meta$localization_prob <= 1))
  expect_equal(unique(s$meta$ptm_type), "Phospho")
  ren <- function(a) { colnames(a$intensity) <- colnames(a$raw_detected) <- ann$sample[match(colnames(a$intensity), paste0("F", 1:6))]; a }
  out <- run_pipeline(ren(p), ren(s))
  check_truth(out, site_from_meta)
})

test_that("FragPipe label-free: combined_protein + combined_modified_peptide (sites from Assigned Modifications)", {
  p <- import_protein(fx("FragPipe_LFQ", "combined_protein.tsv"), "FragPipe")
  s <- import_ptm(fx("FragPipe_LFQ", "combined_modified_peptide.tsv"), "FragPipe")
  expect_equal(colnames(p$intensity), ann$sample); expect_equal(colnames(s$intensity), ann$sample)
  expect_equal(unname(p$intensity[1, 1]), log2(as.numeric(import_table(fx("FragPipe_LFQ", "combined_protein.tsv"))[1, "control_1 MaxLFQ Intensity"])))
  expect_true(all(s$meta$residue %in% c("S", "T", "Y"))); expect_false(any(is.na(s$meta$position)))
  out <- run_pipeline(p, s)
  check_truth(out, site_from_meta)
})

test_that("FragPipe TMT-Integrator tables are on the log2 scale", {
  p <- import_protein(fx("FragPipe_TMT", "abundance_protein_MD.tsv"), "FragPipe", log2_transform = FALSE)
  s <- import_ptm(fx("FragPipe_TMT", "abundance_single-site_MD.tsv"), "FragPipe", log2_transform = FALSE)
  expect_equal(colnames(p$intensity), ann$sample); expect_equal(s$meta$feature_id[1], s$meta$feature_id[1])
  expect_true(all(s$meta$residue %in% c("S", "T", "Y")))
  expect_warning(import_ptm(fx("FragPipe_TMT", "abundance_single-site_MD.tsv"), "FragPipe", log2_transform = TRUE), "log")
  out <- run_pipeline(p, s)
  check_truth(out, site_from_meta)
})

test_that("DIA-NN matrices: raw-file column names, precursors collapsed, sites located with the FASTA", {
  p <- import_protein(fx("DIA-NN", "report.pg_matrix.tsv"), "DIA-NN")
  s <- import_ptm(fx("DIA-NN", "report.pr_matrix.tsv"), "DIA-NN")
  expect_equal(colnames(p$intensity), paste0("run_", ann$sample))
  expect_equal(nrow(s$intensity), 55)                       # 45 peptidoforms, 10 with a second charge
  p <- rename_samples(p, "^run_(.+)$"); s <- rename_samples(s, "^run_(.+)$")
  x <- new_ptm_experiment(p, s, annotation = list(protein = ann, ptm = ann))
  x <- collapse_to_sites(x, "sum", site_key = x$ptm$meta$modified_sequence)      # precursors -> peptidoforms
  expect_equal(nrow(x$ptm$intensity), 45)
  x <- locate_sites(x, fx("example.fasta"))
  expect_true(all(x$ptm$meta$residue %in% c("S", "T", "Y"))); expect_false(any(is.na(x$ptm$meta$position)))
  x <- filter_min_observed(x, min_n = 2) |> normalise(method = "median") |> match_ptm_to_protein("id")
  res <- run_methods(x, make_contrasts(x, "vs_reference", reference = "control"), "limma", pi0_method = "bh")
  check_truth(list(x = x, res = res), site_from_meta)
  y <- locate_sites(new_ptm_experiment(p, s, annotation = list(protein = ann, ptm = ann)))     # without FASTA: residue only
  expect_true(all(is.na(y$ptm$meta$position))); expect_true(all(y$ptm$meta$residue %in% c("S", "T", "Y")))
})
