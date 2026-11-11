#' UpSet plot of the intersections of sets
#'
#' Bars count the elements that belong to exactly the sets marked by the dots below them (an
#' UpSet plot, which reads like a Venn diagram for any number of sets); the size of each set is
#' given next to its name.
#' @param S logical matrix, one row per element and one column per set (column names = set names).
#' @param title,subtitle plot title and subtitle.
#' @param max_bars largest number of intersections shown (the largest ones).
#' @return a ggplot object.
#' @export
plot_upset <- function(S, title = NULL, subtitle = NULL, max_bars = 20) {
  S <- as.matrix(S); storage.mode(S) <- "logical"; S[is.na(S)] <- FALSE
  sets <- colnames(S); k <- length(sets)
  S <- S[rowSums(S) > 0, , drop = FALSE]
  if (!nrow(S)) return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = "Nothing to show") +
                         ggplot2::theme_void() + ggplot2::labs(title = title))
  key <- apply(S, 1, function(v) paste(as.integer(v), collapse = ""))
  tab <- sort(table(key), decreasing = TRUE)
  shown <- utils::head(tab, max_bars)
  if (length(tab) > max_bars) subtitle <- paste0(if (!is.null(subtitle)) paste0(subtitle, "\n") else "",
                                                 sprintf("the %d largest of %d intersections", max_bars, length(tab)))
  combos <- names(shown); n <- as.integer(shown); top <- max(n); step <- top * if (k > 8) 0.07 else 0.15
  x <- seq_along(combos)
  dots <- do.call(rbind, lapply(x, function(i) data.frame(x = i, y = -step * seq_len(k), on = strsplit(combos[i], "")[[1]] == "1")))
  segs <- do.call(rbind, lapply(x, function(i) { d <- dots[dots$x == i & dots$on, ]
    if (nrow(d) < 2) NULL else data.frame(x = i, y = min(d$y), yend = max(d$y)) }))
  lab <- sprintf("%s (%d)", sets, colSums(S))
  br <- pretty(c(0, top)); br <- br[br <= top]
  fs <- if (k > 20) 6 else if (k > 10) 8 else 10
  p <- ggplot2::ggplot() +
    ggplot2::geom_col(data = data.frame(x = x, n = n), ggplot2::aes(x = .data$x, y = .data$n), fill = "#3b6ea5", width = 0.6) +
    ggplot2::geom_text(data = data.frame(x = x, n = n), ggplot2::aes(x = .data$x, y = .data$n, label = .data$n), vjust = -0.4,
                       size = if (length(x) > 12) 2.6 else 3.4) +
    ggplot2::geom_point(data = dots, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$on), size = if (k > 10) 1.8 else 3.2, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#222222", `FALSE` = "#dcdcdc"))
  if (!is.null(segs)) p <- p + ggplot2::geom_segment(data = segs, ggplot2::aes(x = .data$x, xend = .data$x, y = .data$y, yend = .data$yend),
                                                     linewidth = if (k > 10) 0.6 else 1, colour = "#222222")
  p + ggplot2::scale_y_continuous(breaks = c(-step * seq_len(k), br), labels = c(lab, br), expand = ggplot2::expansion(mult = c(0.03, 0.1))) +
    ggplot2::scale_x_continuous(breaks = NULL, limits = c(0.4, length(x) + 0.6)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#888888", linewidth = 0.3) +
    ggplot2::labs(x = NULL, y = NULL, title = title, subtitle = subtitle) + theme_diffptm(11) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = fs, colour = "#222222"))
}

#' Venn diagram of two to four sets
#' @inheritParams plot_upset
#' @return a ggplot object (an UpSet plot when there are more than four sets).
#' @export
plot_venn <- function(S, title = NULL, subtitle = NULL) {
  S <- as.matrix(S); storage.mode(S) <- "logical"; S[is.na(S)] <- FALSE
  k <- ncol(S); if (k < 2 || k > 4) return(plot_upset(S, title, subtitle))
  lay <- switch(as.character(k),
    `2` = list(c(-0.55, 0, 1, 1, 0), c(0.55, 0, 1, 1, 0)),
    `3` = list(c(-0.55, 0.35, 1, 1, 0), c(0.55, 0.35, 1, 1, 0), c(0, -0.6, 1, 1, 0)),
    `4` = list(c(-0.554, 0.191, 1.376, 0.634, 0.739), c(-0.084, -0.044, 1.376, 0.634, 0.739),
               c(0.084, -0.044, 1.376, 0.634, pi - 0.739), c(0.554, 0.191, 1.376, 0.634, pi - 0.739)))
  inell <- function(px, py, e) { dx <- px - e[1]; dy <- py - e[2]; u <- dx * cos(e[5]) + dy * sin(e[5]); v <- -dx * sin(e[5]) + dy * cos(e[5]); (u / e[3])^2 + (v / e[4])^2 <= 1 }
  outline <- do.call(rbind, lapply(seq_len(k), function(i) { e <- lay[[i]]; t <- seq(0, 2 * pi, length.out = 240)
    data.frame(set = colnames(S)[i], x = e[1] + e[3] * cos(t) * cos(e[5]) - e[4] * sin(t) * sin(e[5]),
               y = e[2] + e[3] * cos(t) * sin(e[5]) + e[4] * sin(t) * cos(e[5])) }))
  # label of each region at its most interior grid point
  g <- expand.grid(x = seq(min(outline$x), max(outline$x), length.out = 140), y = seq(min(outline$y), max(outline$y), length.out = 140))
  M <- sapply(lay, function(e) inell(g$x, g$y, e)); rk <- apply(M, 1, function(v) paste(as.integer(v), collapse = ""))
  cnt <- table(factor(apply(S, 1, function(v) paste(as.integer(v), collapse = "")), levels = unique(rk)))
  labs <- do.call(rbind, lapply(setdiff(unique(rk), strrep("0", k)), function(r) {
    inside <- which(rk == r); out <- which(rk != r)
    ob <- out[sample.int(length(out), min(length(out), 3000))]
    d <- vapply(inside, function(i) min((g$x[ob] - g$x[i])^2 + (g$y[ob] - g$y[i])^2), 0)
    j <- inside[which.max(d)]; data.frame(x = g$x[j], y = g$y[j], n = as.integer(cnt[[r]]))
  }))
  cen <- do.call(rbind, lapply(seq_len(k), function(i) { o <- outline[outline$set == colnames(S)[i], ]
    j <- which.max(sqrt((o$x - mean(o$x))^2 + (o$y - mean(o$y))^2) * ifelse(o$y >= mean(o$y), 1, 0.999))
    data.frame(set = colnames(S)[i], x = o$x[j], y = o$y[j]) }))
  pal <- palette_diffptm(k)
  ggplot2::ggplot() +
    ggplot2::geom_polygon(data = outline, ggplot2::aes(x = .data$x, y = .data$y, group = .data$set, fill = .data$set), alpha = 0.18, colour = NA) +
    ggplot2::geom_path(data = outline, ggplot2::aes(x = .data$x, y = .data$y, group = .data$set, colour = .data$set), linewidth = 0.6) +
    ggplot2::geom_text(data = labs, ggplot2::aes(x = .data$x, y = .data$y, label = .data$n), size = if (k == 4) 3.2 else 4) +
    ggplot2::scale_fill_manual(values = pal, name = NULL, labels = sprintf("%s (%d)", colnames(S), colSums(S))) +
    ggplot2::scale_colour_manual(values = pal, name = NULL, labels = sprintf("%s (%d)", colnames(S), colSums(S))) +
    ggplot2::coord_equal() + ggplot2::labs(title = title, subtitle = subtitle) + theme_diffptm(11) +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(), legend.position = "bottom") +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2), colour = ggplot2::guide_legend(ncol = 2))
}

#' Detection of features in each sample, as a logical matrix
#' @param x a `ptm_experiment`.
#' @param assay `"protein"` or `"ptm"`.
#' @param level `"feature"`, or `"protein"` to count, for the PTM assay, the proteins carrying at
#'   least one detected modified feature (first accession of each feature).
#' @return logical matrix (features or proteins x samples), detection before imputation.
#' @export
detection_sets <- function(x, assay = "ptm", level = c("feature", "protein")) {
  level <- match.arg(level); D <- x[[assay]]$raw_detected
  if (level == "protein") {
    acc <- vapply(split_ids(x[[assay]]$meta$protein_ids, x[[assay]]$sep, strip_isoforms = FALSE), function(v) v[1] %||% NA_character_, "")
    acc[is.na(acc)] <- x[[assay]]$meta$feature_id[is.na(acc)]
    D <- rowsum(D * 1, acc) > 0
  }
  D
}

#' Overlap of the features detected in the samples
#'
#' An UpSet plot for all samples, or for the samples of one condition; a Venn diagram when the
#' condition has two to four samples. Features absent from every sample shown are left out.
#' @inheritParams detection_sets
#' @param group a group (condition) whose samples are compared; `NULL` = every sample.
#' @param type `"auto"` (Venn for 2-4 samples, UpSet otherwise), `"venn"` or `"upset"`.
#' @return a ggplot object.
#' @export
plot_detection_overlap <- function(x, assay = "ptm", group = NULL, level = c("feature", "protein"), type = c("auto", "venn", "upset")) {
  level <- match.arg(level); type <- match.arg(type)
  D <- detection_sets(x, assay, level)
  what <- if (assay == "protein") "proteins" else if (level == "protein") "modified proteins" else "PTM features"
  if (!is.null(group)) {
    g <- as.character(groups_of(x, assay)); D <- D[, g == group, drop = FALSE]
  }
  title <- sprintf("Detected %s%s", what, if (is.null(group)) ", all samples" else paste(", condition", group))
  sub <- sprintf("%d %s detected in at least one of these samples, %d in all", sum(rowSums(D) > 0), what, sum(rowSums(D) == ncol(D)))
  if (type == "venn" || (type == "auto" && ncol(D) >= 2 && ncol(D) <= 4)) plot_venn(D, title, sub) else plot_upset(D, title, sub)
}

#' Hierarchical clustering of the samples
#'
#' Ward's method (`ward.D2`) on 1 - Pearson correlation of the log2 intensities, or on 1 - Jaccard
#' index of the detection patterns; leaves are coloured by condition.
#' @inheritParams plot_pca
#' @param measure `"pearson"` or `"jaccard"`.
#' @export
plot_dendrogram <- function(x, assay = "protein", measure = c("pearson", "jaccard")) {
  measure <- match.arg(measure)
  sim <- if (measure == "pearson") qc_pearson(x, assay) else qc_jaccard(x, assay)
  h <- qc_dendrogram(sim)
  n <- length(h$order); xpos <- numeric(n); xpos[h$order] <- seq_len(n)
  node_x <- numeric(nrow(h$merge)); node_y <- h$height
  pos <- function(i) if (i < 0) xpos[-i] else node_x[i]
  hy <- function(i) if (i < 0) 0 else node_y[i]
  seg <- vector("list", nrow(h$merge))
  for (i in seq_len(nrow(h$merge))) {
    a <- h$merge[i, 1]; b <- h$merge[i, 2]; xa <- pos(a); xb <- pos(b); node_x[i] <- (xa + xb) / 2
    seg[[i]] <- data.frame(x = c(xa, xb, xa), xend = c(xa, xb, xb), y = c(hy(a), hy(b), node_y[i]), yend = c(node_y[i], node_y[i], node_y[i]))
  }
  seg <- do.call(rbind, seg)
  smp <- rownames(sim); g <- stats::setNames(as.character(groups_of(x, assay)), x$annotation[[assay]]$sample)
  leaves <- data.frame(x = xpos, sample = smp, group = g[smp])
  fs <- label_size(n) * 0.36
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg, ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend), linewidth = 0.4, colour = "#444444") +
    ggplot2::geom_point(data = leaves, ggplot2::aes(x = .data$x, y = 0, colour = .data$group), size = if (n > 40) 1.2 else 2.2) +
    ggplot2::geom_text(data = leaves, ggplot2::aes(x = .data$x, y = -max(seg$yend) * 0.02, label = .data$sample, colour = .data$group),
                       hjust = 1, size = fs, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = palette_diffptm(length(unique(leaves$group))), name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.35, 0.05)), breaks = function(l) { b <- pretty(c(0, l[2])); b[b >= 0] }) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = sprintf("height (Ward's method, distance 1 - %s)", if (measure == "pearson") "Pearson correlation" else "Jaccard index"),
                  title = sprintf("Clustering of the samples, %s (%s)", if (assay == "protein") "total proteome" else "PTM",
                                  if (measure == "pearson") "intensities" else "detection")) +
    theme_diffptm(10) + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                                       panel.grid.major.y = ggplot2::element_blank(), legend.position = "bottom")
}

#' Variance explained by the principal components (scree plot)
#' @inheritParams plot_pca
#' @param n_pc number of components shown.
#' @export
plot_pca_scree <- function(x, assay = "protein", n_pc = 10) {
  M <- x[[assay]]$intensity; M <- M[stats::complete.cases(M), , drop = FALSE]
  if (nrow(M) < 3 || ncol(M) < 3) stop("Too few complete features or samples for a PCA")
  pc <- stats::prcomp(t(M), center = TRUE, scale. = FALSE)
  v <- pc$sdev^2 / sum(pc$sdev^2); k <- seq_len(min(n_pc, length(v)))
  d <- data.frame(pc = factor(paste0("PC", k), levels = paste0("PC", k)), v = v[k], cum = cumsum(v)[k])
  ggplot2::ggplot(d, ggplot2::aes(x = .data$pc)) +
    ggplot2::geom_col(ggplot2::aes(y = .data$v), fill = "#3b6ea5", width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(y = .data$v, label = sprintf("%.0f%%", 100 * .data$v)), vjust = -0.4, size = 3) +
    ggplot2::geom_line(ggplot2::aes(y = .data$cum, group = 1), colour = "#b84a3a") + ggplot2::geom_point(ggplot2::aes(y = .data$cum), colour = "#b84a3a", size = 1.6) +
    ggplot2::scale_y_continuous(labels = function(b) paste0(100 * b, "%"), limits = c(0, 1.05)) +
    ggplot2::labs(x = NULL, y = "variance explained", title = sprintf("PCA, variance explained, %s", if (assay == "protein") "total proteome" else "PTM"),
                  subtitle = sprintf("bars: each component; line: cumulative (%d complete features)", nrow(M))) + theme_diffptm(11)
}

#' Mean intensity in the two groups of a contrast, per feature
#'
#' Each point is a feature; points off the diagonal change between the groups; significant
#' features are highlighted. For the usage analysis the PTM-vs-protein fold change plot
#' ([plot_fc_comparison()]) is the corresponding view.
#' @param x a processed `ptm_experiment`.
#' @param res results of one analysis (`protein_DPA` or `ptm_DPA`), one contrast and one method
#'   (columns `feature_id`, `log2FC`, `adj_pvalue`).
#' @param group1,group2 the two groups of the contrast (first minus second).
#' @param assay `"protein"` or `"ptm"`.
#' @param fdr,min_log2fc significance thresholds.
#' @param title plot title.
#' @export
plot_condition_means <- function(x, res, group1, group2, assay = "ptm", fdr = 0.05, min_log2fc = 0, title = NULL) {
  g <- as.character(groups_of(x, assay)); M <- x[[assay]]$intensity
  mean_of <- function(grp) { cols <- unlist(strsplit(grp, "+", fixed = TRUE)); rowMeans(M[, g %in% cols, drop = FALSE], na.rm = TRUE) }
  d <- data.frame(feature_id = rownames(M), m1 = mean_of(group1), m2 = mean_of(group2))
  d$sig <- d$feature_id %in% res$feature_id[!is.na(res$adj_pvalue) & res$adj_pvalue < fdr & abs(res$log2FC) >= min_log2fc]
  d <- d[is.finite(d$m1) & is.finite(d$m2), ]; d <- d[order(d$sig), ]
  ggplot2::ggplot(d, ggplot2::aes(x = .data$m2, y = .data$m1, colour = .data$sig)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "#999999") +
    ggplot2::geom_point(size = 0.9, alpha = 0.7) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "#c8c8c8", `TRUE` = "#2a6fdb"), labels = c("Not significant", "Significant"), name = NULL) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = sprintf("mean log2 intensity, %s", group2), y = sprintf("mean log2 intensity, %s", group1), title = title,
                  subtitle = sprintf("%d of %d features significant (FDR %s)", sum(d$sig), nrow(d), fdr)) +
    theme_diffptm(11) + ggplot2::theme(legend.position = "bottom")
}
