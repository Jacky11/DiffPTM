#' msqrob2PTM method ("correct then model")
#'
#' Adapter following the msqrob2PTM workflow (Demeulemeester et al., MCP 2024):
#' 1. `ptm_DPA`: robust regression with empirical Bayes variance moderation (`msqrob2::msqrob`, `robust = TRUE`, no ridge penalty) on the PTM assay with
#'    formula `~ 0 + group` (+ `subject` when present).
#' 2. `ptm_DPU`: PTM intensities are first corrected by subtracting the matched parent
#'    protein log2 intensity *in the same sample* (requires the PTM and proteome assays to
#'    share sample names — msqrob2PTM's within-sample normalisation), then modelled as in 1.
#'    When sample names differ between assays, the per-group mean protein intensity is
#'    subtracted instead and a warning is issued (this approximation must be reported).
#' 3. `protein_DPA`: msqrob on the protein assay.
#'
#' Contrasts are evaluated with `msqrob2::hypothesisTest`.
#'
#' @inheritParams run_methods
#' @param robust passed to `msqrob2::msqrob()`: M-estimation (default) or ordinary least squares.
#' @param ... ignored.
#' @export
method_msqrob2ptm <- function(x, contrasts, analyses = c("protein_DPA", "ptm_DPA", "ptm_DPU"), robust = TRUE, ...) {
  for (p in c("msqrob2", "QFeatures", "SummarizedExperiment"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Package ", p, " is not installed")

  fit_one <- function(M, ann, analysis) {
    g <- factor(if (!all(is.na(ann$factor))) paste(ann$condition, ann$factor, sep = ".") else ann$condition)
    cd <- data.frame(group = g, row.names = colnames(M))
    use_subject <- !all(is.na(ann$subject))
    if (use_subject) cd$subject <- factor(ann$subject)
    use_mix <- !use_subject && !is.null(ann$mixture) && length(unique(stats::na.omit(ann$mixture))) > 1
    if (use_mix) cd$mixture <- factor(ann$mixture)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(intensity = M), colData = cd)
    qf <- QFeatures::QFeatures(list(features = se), colData = cd)
    form <- if (use_subject) ~ 0 + group + subject else if (use_mix) ~ 0 + group + mixture else ~ 0 + group
    qf <- msqrob2::msqrob(qf, i = "features", formula = form, robust = robust)
    lv <- stats::setNames(paste0("group", levels(g)), levels(g))
    gexpr <- function(s) {
      m <- trimws(strsplit(s, "+", fixed = TRUE)[[1]])
      if (length(m) == 1) lv[m] else paste0("(", paste(lv[m], collapse = " + "), ")/", length(m))
    }
    exprs <- mapply(function(a, b) paste0(gexpr(a), " - ", gexpr(b), " = 0"), contrasts$group1, contrasts$group2, USE.NAMES = FALSE)
    L <- msqrob2::makeContrast(exprs, parameterNames = unname(lv))
    qf <- msqrob2::hypothesisTest(qf, i = "features", contrast = L)
    rd <- SummarizedExperiment::rowData(qf[["features"]])
    out <- lapply(seq_len(ncol(L)), function(k) {
      r <- rd[[colnames(L)[k]]]
      data.frame(feature_id = rownames(M), analysis = analysis, contrast = contrasts$name[k],
                 log2FC = r$logFC, se = r$se, df = r$df, statistic = r$t, pvalue = r$pval,
                 adj_pvalue = r$adjPval, method = "msqrob2ptm", stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
  }

  res <- list()
  if ("protein_DPA" %in% analyses && !is.null(x$protein))
    res$prot <- fit_one(x$protein$intensity, x$annotation$protein, "protein_DPA")
  if ("ptm_DPA" %in% analyses && !is.null(x$ptm))
    res$dpa <- fit_one(x$ptm$intensity, x$annotation$ptm, "ptm_DPA")
  if ("ptm_DPU" %in% analyses && !is.null(x$ptm)) {
    if (is.null(x$ptm$protein_intensity)) stop("Run match_ptm_to_protein() first")
    P <- x$ptm$protein_intensity
    common <- intersect(colnames(x$ptm$intensity), colnames(P))
    if (length(common) == ncol(x$ptm$intensity)) {
      corrected <- x$ptm$intensity - P[, colnames(x$ptm$intensity), drop = FALSE]
    } else {
      warning("PTM and proteome samples differ: protein correction uses per-group mean protein intensity")
      gp <- groups_of(x, "protein"); gq <- as.character(groups_of(x, "ptm"))
      Pg <- sapply(levels(gp), function(lv) rowMeans(P[, gp == lv, drop = FALSE], na.rm = TRUE))
      corrected <- x$ptm$intensity - Pg[, gq, drop = FALSE]
    }
    res$dpu <- fit_one(corrected, x$annotation$ptm, "ptm_DPU")
  }
  out <- do.call(rbind, res)
  if (is.null(out)) return(empty_result())
  rownames(out) <- NULL
  out
}
