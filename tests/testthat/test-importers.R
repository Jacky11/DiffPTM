write_tsv <- function(df) { f <- tempfile(fileext = ".txt"); data.table::fwrite(df, f, sep = "\t", na = ""); f }

test_that("MaxQuant proteinGroups: LFQ preferred over Intensity, flags parsed, zeros -> NA", {
  d <- data.frame(check.names = FALSE,
    `Majority protein IDs` = c("P1", "P2;P3", "REV__P4", "CON__P5"), `Gene names` = c("g1", "g2", "", "g5"),
    `Intensity` = 1:4 * 1e6, `Intensity A1` = c(1e6, 2e6, 3e6, 0), `Intensity B1` = c(1e6, 2e6, 3e6, 4e6),
    `LFQ intensity A1` = c(5e6, 0, 7e6, 8e6), `LFQ intensity B1` = c(5e6, 6e6, 7e6, 8e6),
    Reverse = c("", "", "+", ""), `Potential contaminant` = c("", "", "", "+"))
  a <- import_protein(write_tsv(d), "MaxQuant")
  expect_equal(colnames(a$intensity), c("A1", "B1"))
  expect_equal(a$intensity[1, "A1"], log2(5e6))
  expect_true(is.na(a$intensity[2, "A1"]))
  expect_equal(a$meta$Reverse, c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(a$meta$Potential.contaminant, c(FALSE, FALSE, FALSE, TRUE))
  x <- new_ptm_experiment(a, NULL)
  expect_equal(nrow(filter_flagged(x, "protein")$protein$intensity), 2)
})

test_that("MaxQuant sites table: multiplicity expanded / first / sum, totals ignored, localization kept", {
  d <- data.frame(check.names = FALSE,
    Proteins = c("P1", "P2"), `Positions within proteins` = c("10", "20;25"), `Amino acid` = c("S", "T"), Position = c(10, 20),
    `Sequence window` = c("AAASAAA", "AAATAAA"), `Localization prob` = c(0.99, 0.5),
    `Intensity A1` = c(300, 30), `Intensity B1` = c(300, 30),
    `Intensity A1___1` = c(100, 10), `Intensity A1___2` = c(200, 0), `Intensity A1___3` = c(0, 20),
    `Intensity B1___1` = c(100, 10), `Intensity B1___2` = c(200, 0), `Intensity B1___3` = c(0, 20),
    Reverse = c("", ""), `Potential contaminant` = c("", ""))
  f <- write_tsv(d)
  a <- import_ptm(f, "MaxQuant", multiplicity = "expand")
  expect_equal(colnames(a$intensity), c("A1", "B1"))
  expect_equal(nrow(a$intensity), 4)                      # P1 x {1,2}, P2 x {1,3}; all-NA rows dropped
  expect_setequal(a$meta$multiplicity, c("1", "1", "2", "3"))
  expect_equal(a$meta$localization_prob[1], 0.99)
  b <- import_ptm(f, "MaxQuant", multiplicity = "first")
  expect_equal(nrow(b$intensity), 2); expect_equal(unname(b$intensity[1, "A1"]), log2(100))
  s <- import_ptm(f, "MaxQuant", multiplicity = "sum")
  expect_equal(unname(s$intensity[1, "A1"]), log2(300))
})

test_that("Spectronaut PTM pivot: per-sample site probability -> max, quant columns parsed", {
  d <- data.frame(check.names = FALSE,
    PTM.ProteinId = c("P1", "P2"), PTM.CollapseKey = c("P1_S10_M1", "P2_T5_M1"), PTM.SiteLocation = c(10, 5),
    PTM.ModificationTitle = "Phospho (STY)", PEP.StrippedSequence = c("AAASK", "LLTK"),
    `[1] run_A1.raw.PTM.SiteProbability` = c(0.9, 0.2), `[2] run_B1.raw.PTM.SiteProbability` = c(0.6, 0.8),
    `[1] run_A1.raw.PTM.Quantity` = c(1000, 2000), `[2] run_B1.raw.PTM.Quantity` = c(1000, 2000))
  a <- import_ptm(write_tsv(d), "Spectronaut")
  expect_equal(colnames(a$intensity), c("run_A1.raw", "run_B1.raw"))
  expect_equal(a$meta$localization_prob, c(0.9, 0.8))
  expect_equal(a$meta$ptm_type[1], "Phospho (STY)")
})

test_that("FragPipe: MaxLFQ preferred, single-site index parsed into residue/position", {
  p <- data.frame(check.names = FALSE, Protein = c("sp|P1|X", "sp|P2|Y"), `Protein ID` = c("P1", "P2"), Gene = c("g1", "g2"),
                  `A1 Spectral Count` = 1:2, `A1 Intensity` = c(10e5, 20e5), `A1 MaxLFQ Intensity` = c(11e5, 21e5),
                  `B1 Spectral Count` = 1:2, `B1 Intensity` = c(10e5, 20e5), `B1 MaxLFQ Intensity` = c(11e5, 21e5))
  a <- import_protein(write_tsv(p), "FragPipe")
  expect_equal(colnames(a$intensity), c("A1", "B1")); expect_equal(unname(a$intensity[1, 1]), log2(11e5))
  expect_warning(import_protein(write_tsv(within(p, { `A1 MaxLFQ Intensity` <- c(11, 21); `B1 MaxLFQ Intensity` <- c(11, 21) })), "FragPipe"), "log2")
  s <- data.frame(check.names = FALSE, Index = c("P1_S10", "P2_K77"), Gene = c("g1", "g2"), ProteinID = c("P1", "P2"),
                  Peptide = c("AAASK", "LLKK"), SequenceWindow = "x", Start = 1, End = 5, MaxPepProb = c(0.99, 0.8),
                  ReferenceIntensity = 1, A1 = c(20, 21), B1 = c(20.5, 21.5))
  b <- import_ptm(write_tsv(s), "FragPipe", log2_transform = FALSE)
  expect_equal(colnames(b$intensity), c("A1", "B1"))
  expect_equal(b$meta$residue, c("S", "K")); expect_equal(b$meta$position, c("10", "77"))
  expect_equal(b$meta$feature_id, c("P1_S10", "P2_K77"))
})

test_that("Generic format round-trips and column override works", {
  d <- data.frame(check.names = FALSE, protein_id = c("P1", "P2"), gene = c("a", "b"), `intensity_A1` = c(1, 2), `intensity_B1` = c(3, 4))
  a <- import_protein(write_tsv(d), "Generic", log2_transform = FALSE)
  expect_equal(unname(a$intensity[2, "B1"]), 4)
  d2 <- data.frame(check.names = FALSE, protein_id = "P1", gene = "a", `abund.X` = 5, `abund.Y` = 6)
  b <- import_protein(write_tsv(d2), "Generic", quant_pattern = "^abund\\.(.+)$", log2_transform = FALSE)
  expect_equal(colnames(b$intensity), c("X", "Y"))
})

test_that("large MaxQuant intensities (> 2^31) are read as doubles, not integer64", {
  d <- data.frame(check.names = FALSE, `Majority protein IDs` = c("P1", "P2"), `Gene names` = c("g1", "g2"),
                  `LFQ intensity A1` = c(1679400000, 86944000), `LFQ intensity B1` = c(2500000000, 0), Reverse = "", `Potential contaminant` = "")
  a <- import_protein(write_tsv(d), "MaxQuant")
  expect_equal(unname(a$intensity[1, "A1"]), log2(1679400000))
  expect_equal(unname(a$intensity[1, "B1"]), log2(2500000000))
  expect_true(is.na(a$intensity[2, "B1"]))
})

test_that("Spectronaut pivot report with semicolons, decimal commas, 'Filtered' and multi-line quoted cells", {
  f <- tempfile(fileext = ".csv")
  writeLines(c(
    'PG.ProteinGroups;PG.Genes;PG.MolecularWeight;[1] Proteome-01_NS_R1.PG.Quantity;[2] Proteome-02_T_R1.PG.Quantity',
    '"A0A024QZ33;Q9H0G5";NSRP1;"60441,39\n66390,05";7,607966966951104E+05;Filtered',
    'P35237;SERPINB6;42,6;123456,5;234567'), f)
  a <- import_protein(f, "Spectronaut")
  expect_equal(colnames(a$intensity), c("Proteome-01_NS_R1", "Proteome-02_T_R1"))
  expect_equal(unname(a$intensity[1, 1]), log2(760796.6966951104))
  expect_true(is.na(a$intensity[1, 2]))
  expect_equal(unname(a$intensity[2, 2]), log2(234567))
  expect_equal(a$meta$protein_ids[1], "A0A024QZ33;Q9H0G5")
  s <- tempfile(fileext = ".csv")
  writeLines(c(
    'PG.Genes;PTM.ProteinId;PTM.CollapseKey;PTM.Multiplicity;PTM.ModificationTitle;PTM.SiteAA;PTM.SiteLocation;PTM.FlankingRegion;[1] Phospho-01_NS_R1.PTM.SiteProbability;[2] Phospho-02_T_R1.PTM.SiteProbability;[1] Phospho-01_NS_R1.PTM.Quantity;[2] Phospho-02_T_R1.PTM.Quantity',
    'NSRP1;A0A024QZ33;A0A024QZ33_T221_M1;1;Phospho (STY);T;221;RREKVIETPENDFKH;0,99;Filtered;1234,5;2345,5',
    'NSRP1;Q9H0G5;Q9H0G5_T275_M1;1;Phospho (STY);T;275;RREKVIETPENDFKH;0,99;Filtered;1234,5;2345,5'), s)
  b <- import_ptm(s, "Spectronaut")
  expect_equal(b$meta$residue, c("T", "T")); expect_equal(b$meta$position, c("221", "275"))
  expect_equal(b$meta$localization_prob, c(0.99, 0.99)); expect_equal(b$meta$multiplicity, c("1", "1"))
  expect_equal(b$meta$label, c("NSRP1 T221", "NSRP1 T275")); expect_equal(b$meta$feature_id, c("A0A024QZ33_T221_M1", "Q9H0G5_T275_M1"))
  expect_equal(unname(b$intensity[1, 1]), log2(1234.5))
  x <- new_ptm_experiment(a, b, annotation = list(protein = make_annotation(colnames(a$intensity), c("NS", "T")), ptm = make_annotation(colnames(b$intensity), c("NS", "T"))))
  expect_equal(nrow(merge_redundant_ptm(x)$ptm$intensity), 1)
})

test_that("rename_samples selects and renames; guess_annotation accepts ranges and -1", {
  m <- matrix(1:12, 3, dimnames = list(NULL, c("mock_1", "Phospho_mock_1", "Phospho_mock_2", "Phospho_virus_24h_1")))
  a <- new_assay(m, data.frame(feature_id = c("a", "b", "c")))
  r <- rename_samples(a, "^Phospho_(.+)$")
  expect_equal(colnames(r$intensity), c("mock_1", "mock_2", "virus_24h_1"))
  expect_error(rename_samples(a, "^(.+)_\\d$"), "duplicated")
  g <- guess_annotation(colnames(r$intensity), "[_]", "-1")
  expect_equal(g$condition, c("mock", "mock", "virus_24h"))
  expect_equal(guess_annotation(c("CXCL12_3min_R1", "NS_R1"), "[_]", "1-2")$condition, c("CXCL12_3min", "NS_R1"))
  expect_equal(keep_samples(a, c("mock_1", "Phospho_mock_2"))$intensity |> colnames(), c("mock_1", "Phospho_mock_2"))
})
