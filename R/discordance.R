#' Detection-profile discordance between a PTM feature and its parent protein
#'
#' Qualitative counterpart of the DPU test (Giai Gianetto, Methods Mol Biol 2023;2426:267-302,
#' section 3.8; the original app's `pepprot_dect*` sheets):
#' within each group, a PTM feature is *discordant* when it is observed in samples where
#' the parent protein is not (`ptm_only`), or the protein is observed while the PTM is not
#' (`protein_only`). Both are computed on the raw detection pattern (before imputation).
#'
#' Requires sample correspondence between assays: by identical sample names, or by
#' `subject` when names differ; otherwise the comparison is done at the group level
#' (counts of observed samples per group).
#'
#' @param x a `ptm_experiment` after [match_ptm_to_protein()].
#' @param min_discordant minimal number of discordant samples in a group to flag a feature.
#' @return data.frame(feature_id, group, n_samples, n_ptm_only, n_protein_only, n_both, flag),
#'   `flag` being one of `concordant`, `ptm_without_protein`, `protein_without_ptm`,
#'   `no_parent_protein` (feature not matched to any proteome row).
#' @export
detection_discordance <- function(x, min_discordant = 2) {
  if (is.null(x$ptm$protein_intensity)) stop("Run match_ptm_to_protein() first")
  Dq <- x$ptm$raw_detected
  Dp <- !is.na(x$ptm$protein_intensity)        # parent protein detected (per protein sample)
  ann_q <- x$annotation$ptm; ann_p <- x$annotation$protein
  gq <- groups_of(x, "ptm"); gp <- groups_of(x, "protein")

  # sample correspondence
  map <- if (all(colnames(Dq) %in% colnames(Dp))) match(colnames(Dq), colnames(Dp)) else
    if (!all(is.na(ann_q$subject)) && !all(is.na(ann_p$subject))) match(paste(gq, ann_q$subject), paste(gp, ann_p$subject)) else NULL

  rows <- lapply(levels(gq), function(lv) {
    cq <- which(gq == lv)
    if (!is.null(map) && !any(is.na(map[cq]))) {
      q <- Dq[, cq, drop = FALSE]; p <- Dp[, map[cq], drop = FALSE]
      n_ptm_only <- rowSums(q & !p); n_prot_only <- rowSums(!q & p); n_both <- rowSums(q & p)
    } else {
      cp <- which(gp == lv)
      nq <- rowSums(Dq[, cq, drop = FALSE]); np <- rowSums(Dp[, cp, drop = FALSE])
      n_ptm_only <- pmax(nq - np, 0); n_prot_only <- pmax(np - nq, 0); n_both <- pmin(nq, np)
    }
    data.frame(feature_id = x$ptm$meta$feature_id, group = lv, n_samples = length(cq),
               n_ptm_only = n_ptm_only, n_protein_only = n_prot_only, n_both = n_both,
               flag = ifelse(x$ptm$meta$n_matched_protein == 0, "no_parent_protein",
                      ifelse(n_ptm_only >= min_discordant, "ptm_without_protein",
                             ifelse(n_prot_only >= min_discordant, "protein_without_ptm", "concordant"))),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}
