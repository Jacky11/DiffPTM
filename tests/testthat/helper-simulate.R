# Small simulated experiment with known truth, reused across tests.
# 2 conditions x 3 replicates; proteins P1..P40; one PTM feature per protein for P1..P30.
# Truth: P1-P5 protein up (log2FC = 2) with PTM following protein (DPA yes, DPU no);
#        P6-P10 PTM up only (log2FC = 2) with flat protein (DPA yes, DPU yes);
#        others null.
simulate_experiment <- function(seed = 1, sd = 0.3, missing = 0) {
  set.seed(seed)
  cond <- rep(c("A", "B"), each = 3)
  samples <- paste0(cond, "_", 1:3)
  n_prot <- 40
  base_p <- stats::rnorm(n_prot, 25, 2)
  P <- matrix(base_p, n_prot, 6) + matrix(stats::rnorm(n_prot * 6, 0, sd), n_prot)
  P[1:5, cond == "B"] <- P[1:5, cond == "B"] + 2
  dimnames(P) <- list(paste0("P", 1:n_prot), samples)
  n_ptm <- 30
  Q <- P[1:n_ptm, ] - 3 + matrix(stats::rnorm(n_ptm * 6, 0, sd), n_ptm)
  Q[6:10, cond == "B"] <- Q[6:10, cond == "B"] + 2
  dimnames(Q) <- list(paste0("P", 1:n_ptm, "_S10"), samples)
  if (missing > 0) { P[sample(length(P), missing * length(P))] <- NA; Q[sample(length(Q), missing * length(Q))] <- NA }
  prot <- new_assay(P, data.frame(feature_id = rownames(P), protein_ids = rownames(P), gene = NA))
  ptm  <- new_assay(Q, data.frame(feature_id = rownames(Q), protein_ids = paste0("P", 1:n_ptm),
                                  stripped_sequence = NA, modified_sequence = NA, residue = "S", position = "10",
                                  localization_prob = 0.99, ptm_type = "Phospho"))
  ann <- make_annotation(samples, cond)
  x <- new_ptm_experiment(prot, ptm, annotation = list(protein = ann, ptm = ann))
  attr(x, "truth") <- list(protein_up = paste0("P", 1:5), ptm_dpa_up = paste0("P", 1:10, "_S10"), ptm_dpu_up = paste0("P", 6:10, "_S10"))
  x
}
