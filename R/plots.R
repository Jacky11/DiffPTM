#' DiffPTM plot theme and palette
#'
#' All figures of the package and of the app share this theme (light grid, no box) and
#' palette: grey for non-significant features, blue for significant ones, one colour per
#' condition for samples.
#' @param base_size base font size.
#' @export
theme_diffptm <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "#e8e8e8", linewidth = 0.4),
                   axis.title = ggplot2::element_text(colour = "#333333"),
                   plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
                   plot.subtitle = ggplot2::element_text(colour = "#555555", size = ggplot2::rel(0.85)),
                   legend.position = "bottom", legend.title = ggplot2::element_text(size = ggplot2::rel(0.9)),
                   strip.text = ggplot2::element_text(face = "bold"))
}

#' @rdname theme_diffptm
#' @param n number of colours.
#' @export
palette_diffptm <- function(n) {
  base <- c("#2a78d6", "#eb6834", "#1baf7a", "#8e5ad6", "#e0a800", "#17a2b8", "#d63384", "#6c757d", "#795548", "#4caf50", "#ff5722", "#607d8b")
  if (n <= length(base)) base[seq_len(n)] else grDevices::hcl.colors(n, "Dark 3")
}

sig_colours <- c("Not significant" = "#c9c9c9", "Significant" = "#2a78d6")

#' Volcano plot for one analysis x contrast
#' @param res tidy results (one method, one analysis, one contrast), optionally with a
#'   `label` column (gene and site) used for the point labels.
#' @param fdr,min_log2fc significance thresholds.
#' @param label column of `res` used for labels of significant features (`"label"` when
#'   present, `feature_id` otherwise).
#' @param max_labels maximal number of labelled points.
#' @param title,subtitle plot title and subtitle.
#' @export
plot_volcano <- function(res, fdr = 0.05, min_log2fc = 0, label = NULL, max_labels = 20, title = NULL, subtitle = NULL) {
  res <- res[!is.na(res$pvalue), ]
  if (is.null(label)) label <- if ("label" %in% names(res)) "label" else "feature_id"
  res$Significant <- ifelse(res$adj_pvalue < fdr & abs(res$log2FC) >= min_log2fc, "Significant", "Not significant")
  res$Significant <- factor(res$Significant, levels = names(sig_colours))
  thr <- suppressWarnings(max(res$pvalue[res$adj_pvalue < fdr], na.rm = TRUE))
  n_sig <- sum(res$Significant == "Significant")
  if (is.null(subtitle)) subtitle <- sprintf("%d of %d features significant (FDR %s%%%s)", n_sig, nrow(res), format(100 * fdr),
                                             if (min_log2fc > 0) sprintf(", |log2FC| >= %s", format(min_log2fc)) else "")
  p <- ggplot2::ggplot(res, ggplot2::aes(x = .data$log2FC, y = -log10(.data$pvalue), colour = .data$Significant)) +
    ggplot2::geom_point(alpha = 0.75, size = 1.4) +
    ggplot2::scale_colour_manual(values = sig_colours, drop = FALSE, name = NULL) +
    ggplot2::labs(x = expression(log[2]~fold~change), y = expression(-log[10]~p), title = title, subtitle = subtitle) +
    theme_diffptm()
  if (min_log2fc > 0) p <- p + ggplot2::geom_vline(xintercept = c(-min_log2fc, min_log2fc), colour = "#888888", linetype = "dashed")
  if (is.finite(thr)) p <- p + ggplot2::geom_hline(yintercept = -log10(thr), colour = "#888888", linetype = "dashed")
  sig <- res[res$Significant == "Significant", ]
  if (nrow(sig) && requireNamespace("ggrepel", quietly = TRUE)) {
    sig <- sig[order(sig$pvalue), ][seq_len(min(max_labels, nrow(sig))), ]
    p <- p + ggrepel::geom_text_repel(data = sig, ggplot2::aes(label = .data[[label]]), size = 2.8, colour = "#222222",
                                      max.overlaps = 20, min.segment.length = 0.2, segment.colour = "#999999", box.padding = 0.4, show.legend = FALSE)
  }
  p
}

#' Histogram of raw p-values
#'
#' A flat histogram with a peak near 0 is what a well-calibrated test on data with some true
#' changes looks like; a histogram rising towards 1 or U-shaped signals a problem.
#' @inheritParams plot_volcano
#' @param res tidy results.
#' @param title plot title.
#' @export
plot_pvalue_histogram <- function(res, title = NULL) {
  res <- res[!is.na(res$pvalue), ]
  ggplot2::ggplot(res, ggplot2::aes(x = .data$pvalue)) +
    ggplot2::geom_histogram(breaks = seq(0, 1, 0.025), fill = "#9cc3ea", colour = "white", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = nrow(res) / 40, colour = "#888888", linetype = "dashed") +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    ggplot2::labs(x = "p-value", y = "features", title = title, subtitle = "dashed line: uniform distribution") + theme_diffptm()
}

#' PTM vs protein log2 fold-change scatter
#'
#' Each point is a PTM feature; its protein change is on the x axis and its own change on the
#' y axis. Points on the diagonal follow their protein (no usage change); points off the
#' diagonal, coloured when significant in the usage (DPU) test, change relative to it.
#'
#' @param res_ptm,res_prot tidy results for `ptm_DPA` and `protein_DPA` on the same
#'   contrast.
#' @param parent named character vector mapping each PTM `feature_id` to its proteome
#'   feature (`setNames(x$ptm$meta$matched_protein, x$ptm$meta$feature_id)`).
#' @param res_dpu optional `ptm_DPU` results used for colouring (significant at `fdr`).
#' @param fdr,min_log2fc thresholds applied to `res_dpu`.
#' @param title plot title.
#' @export
plot_fc_comparison <- function(res_ptm, res_prot, parent, res_dpu = NULL, fdr = 0.05, min_log2fc = 0, title = NULL) {
  first_parent <- sub(";.*$", "", parent[match(res_ptm$feature_id, names(parent))])
  df <- data.frame(feature_id = res_ptm$feature_id, ptm = res_ptm$log2FC,
                   protein = res_prot$log2FC[match(first_parent, res_prot$feature_id)])
  df$label <- if ("label" %in% names(res_ptm)) res_ptm$label else res_ptm$feature_id
  df$usage <- rep("Not significant", nrow(df))
  if (!is.null(res_dpu) && nrow(res_dpu) && nrow(df)) {
    k <- match(df$feature_id, res_dpu$feature_id)
    df$usage <- ifelse(!is.na(res_dpu$adj_pvalue[k]) & res_dpu$adj_pvalue[k] < fdr & abs(res_dpu$log2FC[k]) >= min_log2fc, "Significant", "Not significant")
  }
  df$usage <- factor(df$usage, levels = names(sig_colours))
  df <- df[stats::complete.cases(df[, c("ptm", "protein")]), ]
  df <- df[order(df$usage), ]
  lim <- if (nrow(df)) max(abs(c(df$ptm, df$protein)), na.rm = TRUE) else 1
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$protein, y = .data$ptm, colour = .data$usage)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, colour = "#bbbbbb") +
    ggplot2::geom_hline(yintercept = 0, colour = "#dddddd") + ggplot2::geom_vline(xintercept = 0, colour = "#dddddd") +
    ggplot2::geom_point(alpha = 0.7, size = 1.4) +
    ggplot2::scale_colour_manual(values = sig_colours, drop = FALSE, name = "Usage (DPU)") +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::labs(x = expression(log[2]~FC~(protein)), y = expression(log[2]~FC~(PTM~feature)), title = title,
                  subtitle = "diagonal: the feature follows its protein") + theme_diffptm()
  sig <- df[df$usage == "Significant", ]
  if (nrow(sig) && requireNamespace("ggrepel", quietly = TRUE)) {
    sig <- sig[order(-abs(sig$ptm - sig$protein)), ][seq_len(min(15, nrow(sig))), ]
    p <- p + ggrepel::geom_text_repel(data = sig, ggplot2::aes(label = .data$label), size = 2.8, colour = "#222222",
                                      max.overlaps = 20, segment.colour = "#999999", show.legend = FALSE)
  }
  p
}

#' Heatmap of a similarity / correlation matrix (ggplot2)
#'
#' Samples are ordered by hierarchical clustering when `cluster = TRUE`; axis labels are
#' hidden beyond 40 samples; `groups` adds a row of coloured squares giving the condition
#' of each sample.
#' @param mat square matrix with dimnames.
#' @param title plot title.
#' @param limits colour scale limits.
#' @param cluster reorder samples by clustering (Ward.D2 on 1 - similarity).
#' @param groups optional named vector (sample -> condition).
#' @export
plot_heatmap_matrix <- function(mat, title = "", limits = c(0, 1), cluster = TRUE, groups = NULL) {
  ord <- if (cluster && nrow(mat) > 2) qc_dendrogram(mat)$order else seq_len(nrow(mat))
  mat <- mat[ord, ord, drop = FALSE]
  n <- nrow(mat); samples <- rownames(mat)
  ylev <- if (is.null(groups)) samples else c("condition", samples)
  df <- data.frame(x = factor(rep(samples, times = n), levels = samples), y = factor(rep(samples, each = n), levels = ylev), value = as.vector(mat))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$value)) +
    ggplot2::scale_fill_gradientn(colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"), limits = limits, name = NULL, na.value = "white") +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL, title = title) + theme_diffptm(10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "right", legend.box = "vertical",
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = if (n > 20) 6 else 8),
                   axis.text.y = ggplot2::element_text(size = if (n > 20) 6 else 8))
  if (n > 40) p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) +
    ggplot2::labs(subtitle = sprintf("%d samples, ordered by clustering", n))
  if (!is.null(groups)) {
    g <- data.frame(x = factor(samples, levels = samples), y = factor("condition", levels = ylev), group = as.character(groups[samples]))
    p <- p + ggplot2::geom_point(data = g, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$group), shape = 15, size = if (n > 40) 1.6 else 3) +
      ggplot2::scale_colour_manual(values = palette_diffptm(length(unique(g$group))), name = "condition") +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = if (length(unique(g$group)) > 12) 2 else 1, override.aes = list(size = 3)))
  }
  p
}

#' Number of features observed per sample
#'
#' Runs that quantified far fewer features than the others (failed injections, empty
#' enrichments) stand out; they should be removed in the design step.
#' @inheritParams plot_pca
#' @export
plot_detected_per_sample <- function(x, assay = "ptm") {
  D <- x[[assay]]$raw_detected
  d <- data.frame(sample = factor(colnames(D), levels = colnames(D)), n = colSums(D), group = as.character(groups_of(x, assay)))
  n <- nrow(d)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$sample, y = .data$n, fill = .data$group)) + ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_hline(yintercept = stats::median(d$n) / 2, colour = "#888888", linetype = "dashed") +
    ggplot2::scale_fill_manual(values = palette_diffptm(length(unique(d$group))), name = NULL) +
    ggplot2::labs(x = NULL, y = "features observed", title = sprintf("Features observed per sample, %s", assay),
                  subtitle = "dashed line: half the median; samples below it deserve a look") + theme_diffptm(11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = if (n > 20) 6 else 8))
  if (n > 40) p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  p
}

#' PCA of samples
#' @param x a `ptm_experiment`.
#' @param assay `"protein"` or `"ptm"`.
#' @param label write sample names next to the points (automatic: when at most 30 samples).
#' @export
plot_pca <- function(x, assay = "protein", label = NULL) {
  d <- qc_pca(x, assay); v <- attr(d, "var_explained")
  if (is.null(label)) label <- nrow(d) <= 30
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$PC1, .data$PC2, colour = .data$group)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#dddddd") + ggplot2::geom_vline(xintercept = 0, colour = "#dddddd") +
    ggplot2::geom_point(size = 2.6, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = palette_diffptm(length(unique(d$group))), name = NULL) +
    ggplot2::labs(x = sprintf("PC1 (%.0f%%)", 100 * v[1]), y = sprintf("PC2 (%.0f%%)", 100 * v[2]),
                  title = sprintf("PCA, %s (%d complete features)", assay, attr(d, "n_features") %||% NA)) +
    theme_diffptm() + ggplot2::guides(colour = ggplot2::guide_legend(ncol = min(5, length(unique(d$group))), byrow = TRUE)) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = if (length(unique(d$group)) > 8) 8 else 10))
  if (label && requireNamespace("ggrepel", quietly = TRUE))
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = .data$sample), size = 2.6, colour = "#444444", show.legend = FALSE, max.overlaps = 20)
  p
}

#' Intensity distributions per sample (boxplots)
#' @inheritParams plot_pca
#' @export
plot_intensity_distributions <- function(x, assay = "protein") {
  d <- qc_intensity_distributions(x, assay)
  d$sample <- factor(d$sample, levels = colnames(x[[assay]]$intensity))
  n <- nlevels(d$sample)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$sample, y = .data$log2_intensity, fill = .data$group)) +
    ggplot2::geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = palette_diffptm(length(unique(d$group))), name = NULL) +
    ggplot2::labs(x = NULL, y = "log2 intensity", title = sprintf("Intensity distributions, %s", assay)) + theme_diffptm(11) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = min(6, length(unique(d$group))), byrow = TRUE)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = if (n > 20) 6 else 8))
  if (n > 40) p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  p
}

#' Bar chart of significant features per contrast and method
#' @param res tidy results.
#' @param fdr,min_log2fc thresholds.
#' @param analysis analysis shown.
#' @export
plot_counts <- function(res, fdr = 0.05, min_log2fc = 0, analysis = "ptm_DPU") {
  r <- res[res$analysis == analysis & !is.na(res$adj_pvalue), ]
  r$sig <- r$adj_pvalue < fdr & abs(r$log2FC) >= min_log2fc
  d <- stats::aggregate(sig ~ contrast + method, r, sum)
  d$contrast <- factor(gsub("_vs_", "\nvs\n", d$contrast, fixed = TRUE), levels = gsub("_vs_", "\nvs\n", unique(r$contrast), fixed = TRUE))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$contrast, y = .data$sig, fill = .data$method)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(0.8), width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = .data$sig), position = ggplot2::position_dodge(0.8), vjust = -0.3, size = 3) +
    ggplot2::scale_fill_manual(values = palette_diffptm(length(unique(d$method))), name = NULL) +
    ggplot2::labs(x = NULL, y = sprintf("significant features (FDR %s%%)", format(100 * fdr)), title = analysis) + theme_diffptm() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = if (nlevels(d$contrast) > 4) 7 else 9))
}
