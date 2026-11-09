#' Detection profiles (1 = observed before imputation, 0 = missing)
#' @param x a `ptm_experiment`.
#' @param assay `"protein"` or `"ptm"`.
#' @export
qc_detection <- function(x, assay) {
  d <- x[[assay]]$raw_detected * 1L
  rownames(d) <- x[[assay]]$meta$feature_id
  d
}

#' Jaccard index between samples on detection patterns
#' @inheritParams qc_detection
#' @export
qc_jaccard <- function(x, assay) {
  D <- x[[assay]]$raw_detected
  n <- ncol(D); J <- matrix(NA_real_, n, n, dimnames = list(colnames(D), colnames(D)))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    inter <- sum(D[, i] & D[, j]); uni <- sum(D[, i] | D[, j])
    J[i, j] <- if (uni == 0) NA_real_ else inter / uni
  }
  J
}

#' Pearson correlation between samples (pairwise complete observations)
#' @inheritParams qc_detection
#' @param use passed to [stats::cor()].
#' @export
qc_pearson <- function(x, assay, use = "pairwise.complete.obs") {
  stats::cor(x[[assay]]$intensity, use = use, method = "pearson")
}

#' Hierarchical clustering (Ward.D2) on 1 - similarity
#' @param sim similarity matrix from [qc_jaccard()] or [qc_pearson()].
#' @return an `hclust` object; use `plot()` or ggdendro for display.
#' @export
qc_dendrogram <- function(sim) {
  sim[is.na(sim)] <- 0
  stats::hclust(stats::as.dist(1 - sim), method = "ward.D2")
}

#' PCA of samples on complete features
#' @inheritParams qc_detection
#' @param n_pc number of components returned.
#' @return data.frame(sample, group, PC1, PC2, ...) with `var_explained` attribute.
#' @export
qc_pca <- function(x, assay, n_pc = 3) {
  M <- x[[assay]]$intensity
  M <- M[stats::complete.cases(M), , drop = FALSE]
  if (nrow(M) < 10) {                  # few complete features (many samples, PTM data): impute the sample median on features >= 70 % observed
    M <- x[[assay]]$intensity; M <- M[rowMeans(!is.na(M)) >= 0.7, , drop = FALSE]
    if (nrow(M) < 3) stop("Not enough features observed in at least 70 % of the samples for a PCA")
    med <- apply(M, 2, stats::median, na.rm = TRUE)
    for (j in seq_len(ncol(M))) M[is.na(M[, j]), j] <- med[j]
  }
  p <- stats::prcomp(t(M), center = TRUE, scale. = FALSE)
  k <- min(n_pc, ncol(p$x))
  out <- data.frame(sample = colnames(M), group = as.character(groups_of(x, assay)), p$x[, seq_len(k), drop = FALSE])
  attr(out, "var_explained") <- (p$sdev^2 / sum(p$sdev^2))[seq_len(k)]
  attr(out, "n_features") <- nrow(M)
  out
}

#' Long table of intensities for distribution plots (ridges / boxplots)
#' @inheritParams qc_detection
#' @export
qc_intensity_distributions <- function(x, assay) {
  M <- x[[assay]]$intensity
  g <- groups_of(x, assay)
  data.frame(sample = rep(colnames(M), each = nrow(M)),
             group = rep(as.character(g), each = nrow(M)),
             log2_intensity = as.vector(M))[!is.na(as.vector(M)), ]
}
