###############################################################################
# overlap_clustering.R
#
# Optional analyses of the samples:
#   1 - Pairwise overlap between samples (Jaccard index, overlap coefficient)
#       computed on protein presence/absence, summarised within and between
#       groups, with permutation tests (betadisper / PERMANOVA) that do not
#       depend on group size.
#   2 - Hierarchical clustering (Ward.D2, 1 - Pearson correlation computed on
#       pairwise complete observations, i.e. the sample dendrogram of the report)
#       with multiscale bootstrap resampling (pvclust): AU and BP p-values.
#
# Dependencies: ggplot2, reshape2, openxlsx, officer, pvclust, vegan
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(pvclust)
  library(vegan)
  library(openxlsx)
  library(officer)
})

#------------------------------------------------------------------------------
# Helpers: sample labels and groups
#------------------------------------------------------------------------------

# "Intensity.Ctrl1" / "LFQ.intensity.Ctrl1" -> "Ctrl1"
clean_sample_names <- function(x) {
  x <- sub("^LFQ\\.intensity\\.", "", x)
  x <- sub("^Intensity\\.", "", x)
  x
}

# Group = sample label without its trailing number ("Ctrl1" -> "Ctrl",
# "T.KO3" -> "T.KO"). A manual comma-separated vector can override it.
make_groups <- function(labels, manual = NULL) {
  if (!is.null(manual) && nzchar(trimws(manual))) {
    g <- trimws(unlist(strsplit(manual, ",")))
    if (length(g) != length(labels))
      stop(sprintf("Manual groups: %d labels given for %d samples.", length(g), length(labels)))
  } else {
    # remove the trailing replicate number and any letter suffix:
    # "Ctrl3" -> "Ctrl", "Ctrl_3B" -> "Ctrl", "T.KO_1A" -> "T.KO"
    g <- sub("[._ -]*[0-9]+[A-Za-z]*$", "", labels)
  }
  g <- factor(g, levels = unique(g))
  if (any(table(g) < 2))
    stop(paste0("Each group needs at least 2 samples. Groups found: ",
                paste(names(table(g)), table(g), sep = "=", collapse = ", "),
                ". Fill in the manual groups field (comma-separated, one per sample)."))
  g
}

group_palette <- function(groups) {
  base <- c("#D9D9D9", "#FFFF00", "#DA70D6", "#F07830", "#B5E3A0", "#4CA832",
            "#6BAED6", "#FDAE6B", "#9E9AC8", "#A1D99B", "#FB6A4A", "#BCBDDC")
  lv <- levels(groups)
  setNames(rep(base, length.out = length(lv)), lv)
}

#------------------------------------------------------------------------------
# 1 - Pairwise overlap indices
#------------------------------------------------------------------------------

# P: logical/0-1 matrix proteins x samples (TRUE = protein quantified)
pairwise_overlap <- function(P) {
  P <- (P * 1)
  P[is.na(P)] <- 0
  inter <- crossprod(P)                           # |A n B|
  n     <- colSums(P)                              # |A|
  union <- outer(n, n, "+") - inter                # |A u B|
  list(
    n_prot       = n,
    intersection = inter,
    union        = union,
    jaccard      = inter / union,                  # |A n B| / |A u B|
    overlap      = inter / outer(n, n, pmin),      # Szymkiewicz-Simpson
    sorensen     = 2 * inter / outer(n, n, "+")    # Dice-Sorensen
  )
}

# Long table of all sample pairs (i < j), flagged within / between groups
pairs_table <- function(ov, groups) {
  s  <- colnames(ov$jaccard)
  ij <- which(upper.tri(ov$jaccard), arr.ind = TRUE)
  g1 <- as.character(groups[ij[, 1]]); g2 <- as.character(groups[ij[, 2]])
  data.frame(
    Sample1 = s[ij[, 1]], Sample2 = s[ij[, 2]],
    Group1 = g1, Group2 = g2,
    Type = ifelse(g1 == g2, "Within", "Between"),
    N1 = ov$n_prot[ij[, 1]], N2 = ov$n_prot[ij[, 2]],
    Shared = ov$intersection[ij], Union = ov$union[ij],
    Jaccard = ov$jaccard[ij], Overlap = ov$overlap[ij], Sorensen = ov$sorensen[ij],
    stringsAsFactors = FALSE
  )
}

# Proportion of proteins shared by ALL samples of a set among proteins seen in
# at least one of them (the "% of common proteins"). It decreases
# mechanically with the number of samples, hence the size-matched version.
core_fraction <- function(P, idx) {
  r <- rowSums(P[, idx, drop = FALSE])
  sum(r == length(idx)) / sum(r > 0)
}

size_matched_core <- function(P, idx, k, max_comb = 1000, seed = 123) {
  if (length(idx) <= k) return(c(mean = core_fraction(P, idx), sd = NA, n_subsets = 1))
  cb <- combn(idx, k)
  if (ncol(cb) > max_comb) { set.seed(seed); cb <- cb[, sample(ncol(cb), max_comb)] }
  v <- apply(cb, 2, function(s) core_fraction(P, s))
  c(mean = mean(v), sd = sd(v), n_subsets = ncol(cb))
}

group_summary <- function(P, pt, groups) {
  k <- min(table(groups))
  do.call(rbind, lapply(levels(groups), function(g) {
    idx <- which(groups == g)
    w   <- pt[pt$Type == "Within" & pt$Group1 == g, ]
    b   <- pt[pt$Type == "Between" & (pt$Group1 == g | pt$Group2 == g), ]
    smc <- size_matched_core(P, idx, k)
    data.frame(
      Group = g, n_samples = length(idx), n_pairs = nrow(w),
      Jaccard_within_median = median(w$Jaccard),
      Jaccard_within_Q1 = quantile(w$Jaccard, .25), Jaccard_within_Q3 = quantile(w$Jaccard, .75),
      Jaccard_within_mean = mean(w$Jaccard), Jaccard_within_min = min(w$Jaccard),
      Jaccard_between_median = median(b$Jaccard),
      Overlap_within_median = median(w$Overlap),
      Core_pct_all_samples = 100 * core_fraction(P, idx),
      Core_pct_size_matched = 100 * smc["mean"], Core_pct_size_matched_sd = 100 * smc["sd"],
      Size_matched_k = k, row.names = NULL
    )
  }))
}

# Permutation tests on the Jaccard distance (1 - J)
#  - betadisper + permutest : H0 = same multivariate dispersion (heterogeneity)
#    in all groups; the pairwise version compares groups two by two.
#  - adonis2 (PERMANOVA)    : H0 = same centroid (composition) in all groups.
jaccard_tests <- function(J, groups, nperm = 9999, seed = 123) {
  D  <- as.dist(1 - J)
  bd <- betadisper(D, groups, type = "centroid", bias.adjust = TRUE)
  set.seed(seed)
  pt <- permutest(bd, permutations = nperm, pairwise = TRUE)
  set.seed(seed)
  ad <- adonis2(D ~ groups, permutations = nperm)
  pw <- pt$pairwise$permuted
  pw_tab <- data.frame(Comparison = names(pw), p_perm = as.numeric(pw),
                       p_perm_BH = p.adjust(as.numeric(pw), "BH"))
  list(
    betadisper = bd,
    dispersion_global = data.frame(
      Test = "betadisper (homogeneity of multivariate dispersion)",
      F = pt$tab[1, "F"], df1 = pt$tab[1, "Df"], df2 = pt$tab[2, "Df"],
      p_perm = pt$tab[1, "Pr(>F)"], n_perm = nperm),
    dispersion_pairwise = pw_tab,
    distance_to_centroid = data.frame(Sample = names(bd$distances), Group = groups,
                                      Dist_to_centroid = bd$distances),
    permanova = data.frame(Test = "PERMANOVA (adonis2)", F = ad$F[1],
                           R2 = ad$R2[1], df = ad$Df[1], p_perm = ad$`Pr(>F)`[1],
                           n_perm = nperm)
  )
}

#------------------------------------------------------------------------------
# 1 - Plots
#------------------------------------------------------------------------------

plot_jaccard_heatmap <- function(M, groups, title = "Pairwise Jaccard index", digits = 2) {
  pal <- group_palette(groups)
  o   <- order(groups)
  M   <- M[o, o]; g <- groups[o]
  df  <- reshape2::melt(M, varnames = c("S1", "S2"), value.name = "value")
  df$S1 <- factor(df$S1, levels = rownames(M)); df$S2 <- factor(df$S2, levels = rownames(M))
  n <- nrow(M)
  ggplot(df, aes(S1, S2, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    geom_text(aes(label = formatC(value, format = "f", digits = digits)),
              size = ifelse(n > 20, 1.8, 2.6),
              colour = ifelse(df$value > 0.75, "white", "black")) +
    scale_fill_gradientn(colours = c("white", "lightblue", "darkblue"),
                         limits = c(min(0.3, min(M)), 1), name = "Index") +
    annotate("rect", xmin = seq_len(n) - 0.5, xmax = seq_len(n) + 0.5,
             ymin = -0.6, ymax = 0.4, fill = pal[as.character(g)]) +
    annotate("rect", ymin = seq_len(n) - 0.5, ymax = seq_len(n) + 0.5,
             xmin = -0.6, xmax = 0.4, fill = pal[as.character(g)]) +
    coord_equal(clip = "off") + labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA))
}

plot_within_group <- function(pt, groups, tests = NULL, index = "Jaccard") {
  pal <- group_palette(groups)
  w <- pt[pt$Type == "Within", ]
  w$Group <- factor(w$Group1, levels = levels(groups))
  b <- pt[pt$Type == "Between", ]
  df <- rbind(data.frame(Group = w$Group, value = w[[index]]),
              data.frame(Group = "Between groups", value = b[[index]]))
  df$Group <- factor(df$Group, levels = c(levels(groups), "Between groups"))
  fills <- c(pal, "Between groups" = "grey80")
  sub <- NULL
  if (!is.null(tests)) {
    sub <- sprintf("Dispersion (betadisper): F = %.2f, p = %.4f   |   PERMANOVA: R2 = %.2f, p = %.4f",
                   tests$dispersion_global$F, tests$dispersion_global$p_perm,
                   tests$permanova$R2, tests$permanova$p_perm)
  }
  ggplot(df, aes(Group, value, fill = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
    geom_jitter(width = 0.15, height = 0, size = 1.3, alpha = 0.7) +
    scale_fill_manual(values = fills, guide = "none") +
    labs(x = NULL, y = paste0("Pairwise ", index, " index"),
         title = paste0("Pairwise ", index, " index between samples of the same group"),
         subtitle = sub) +
    theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size = 8))
}

plot_dispersion <- function(tests, groups) {
  pal <- group_palette(groups)
  df <- tests$distance_to_centroid
  df$Group <- factor(df$Group, levels = levels(groups))
  ggplot(df, aes(Group, Dist_to_centroid, fill = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
    geom_jitter(width = 0.12, height = 0, size = 1.5) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(x = NULL, y = "Distance to group centroid (Jaccard distance, PCoA)",
         title = "Within-group heterogeneity (multivariate dispersion)",
         subtitle = sprintf("Permutation test (%d perm.): F = %.2f, p = %.4f",
                            tests$dispersion_global$n_perm, tests$dispersion_global$F,
                            tests$dispersion_global$p_perm)) +
    theme_bw(base_size = 11)
}

#------------------------------------------------------------------------------
# 2 - pvclust with the same distance / linkage as the sample dendrogram of the report
#------------------------------------------------------------------------------

# x: proteins x samples matrix of log2 intensities (NA allowed)
run_pvclust <- function(x, nboot = 10000, method.hclust = "ward.D2",
                        parallel = FALSE, seed = 123) {
  x <- as.matrix(x)
  pvclust(x, method.hclust = method.hclust, method.dist = "correlation",
          use.cor = "pairwise.complete.obs", nboot = nboot,
          parallel = parallel, iseed = seed, quiet = TRUE)
}

# Table of all internal nodes (edges) with AU / BP and their members
pvclust_edge_table <- function(pv, groups) {
  hc <- pv$hclust
  lab <- hc$labels
  members <- vector("list", nrow(hc$merge))
  for (i in seq_len(nrow(hc$merge))) {
    m <- hc$merge[i, ]
    members[[i]] <- unlist(lapply(m, function(k) if (k < 0) lab[-k] else members[[k]]))
  }
  comp <- sapply(members, function(m) {
    tb <- table(droplevels(groups[match(m, lab)]))
    paste(paste0(names(tb), ":", tb), collapse = ", ")
  })
  e <- pv$edges
  data.frame(Edge = seq_len(nrow(e)), Height = hc$height,
             AU_pct = round(100 * e$au), BP_pct = round(100 * e$bp),
             SE_AU = round(e$se.au, 4), n_members = lengths(members),
             Group_composition = comp,
             Members = sapply(members, paste, collapse = ", "),
             Significant_AU95 = e$au >= 0.95)
}

# Horizontal dendrogram with AU (red) and BP
# (green) values at each node; dashed boxes = clusters with AU >= alpha.
plot_pvclust_gg <- function(pv, groups, alpha = 0.95, title = NULL,
                            rescale_height = FALSE) {
  hc  <- pv$hclust
  lab <- hc$labels
  pal <- group_palette(groups)
  n   <- length(lab)
  h   <- hc$height
  if (rescale_height) h <- h / max(h)
  ypos_leaf <- setNames(match(seq_len(n), hc$order), lab)   # leaf -> position
  node_y <- node_x <- numeric(nrow(hc$merge))
  child_xy <- function(k) if (k < 0) c(0, ypos_leaf[-k]) else c(node_x[k], node_y[k])
  seg <- list()
  for (i in seq_len(nrow(hc$merge))) {
    a <- child_xy(hc$merge[i, 1]); b <- child_xy(hc$merge[i, 2])
    node_x[i] <- h[i]; node_y[i] <- (a[2] + b[2]) / 2
    seg[[i]] <- data.frame(x = c(a[1], b[1], h[i]), xend = c(h[i], h[i], h[i]),
                           y = c(a[2], b[2], a[2]), yend = c(a[2], b[2], b[2]))
  }
  seg <- do.call(rbind, seg)
  nodes <- data.frame(x = node_x, y = node_y,
                      AU = round(100 * pv$edges$au), BP = round(100 * pv$edges$bp))
  leaves <- data.frame(label = lab, y = ypos_leaf[lab], group = groups)
  # clusters with AU >= alpha (largest non-overlapping ones, as pvrect)
  pk <- pvpick(pv, alpha = alpha, pv = "au", type = "geq", max.only = TRUE)
  rects <- do.call(rbind, lapply(pk$edges, function(e) {
    yy <- ypos_leaf[pk$clusters[[match(e, pk$edges)]]]
    data.frame(xmin = -0.01 * max(h), xmax = h[e] * 1.03 + 0.01 * max(h),
               ymin = min(yy) - 0.42, ymax = max(yy) + 0.42)
  }))
  dx <- 0.012 * max(h)
  p <- ggplot() +
    geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4)
  if (!is.null(rects))
    p <- p + geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = NA, colour = "red3", linetype = "dashed", linewidth = 0.4)
  p <- p +
    geom_text(data = nodes, aes(x = x - dx, y = y + 0.22, label = AU),
              colour = "red3", size = 2.3, hjust = 1, vjust = 0) +
    geom_text(data = nodes, aes(x = x + dx, y = y + 0.22, label = BP),
              colour = "green4", size = 2.3, hjust = 0, vjust = 0) +
    geom_label(data = leaves, aes(x = -0.02 * max(h), y = y, label = label, fill = group),
               hjust = 1, size = 3, label.size = 0, label.r = unit(0, "pt")) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.18, 0.03))) +
    scale_y_continuous(expand = expansion(add = 0.8)) +
    labs(x = "Height", y = NULL,
         title = title,
         subtitle = sprintf("Numbers at nodes: AU p-value (red) / BP (green), %%, nboot = %d; dashed boxes: AU >= %d%%",
                            pv$nboot[1], round(100 * alpha))) +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 8))
  p
}

#------------------------------------------------------------------------------
# Wrapper called from run_pairwise() (R/pipeline.R)
#------------------------------------------------------------------------------
# x          : proteins x samples matrix of log2 intensities, NA = not quantified
# groups     : factor of sample groups (same order as columns of x)
# prefix     : prefix of output files (e.g. "AllProteins" or "ImmuneMediators")
# doc        : officer pptx object (slides are appended); NULL to skip
# Returns list(doc, files, results)
run_overlap_analyses <- function(x, groups, prefix = "Overlap", doc = NULL,
                                  nboot = 10000, nperm = 9999, parallel = FALSE,
                                  alpha = 0.95, seed = 123, outdir = ".",
                                  modif = FALSE, modif_file = "doc.pptx") {
  x <- as.matrix(x)
  colnames(x) <- clean_sample_names(colnames(x))
  names(groups) <- colnames(x)
  files <- character(0)
  out <- function(f) {
    fn <- paste0(prefix, "_", f)
    if (outdir != ".") fn <- file.path(outdir, fn)
    fn
  }
  save_plot <- function(p, f, w, h) {
    ggsave(out(paste0(f, ".pdf")), p, width = w, height = h, device = cairo_pdf, bg = "white")
    ggsave(out(paste0(f, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
    files <<- c(files, out(paste0(f, ".pdf")), out(paste0(f, ".png")))
    to_doc(p, w, h)
    out(paste0(f, ".png"))
  }
  # editable (vector) copy of each plot appended to doc.pptx, as the other
  # plots of the app when "PPTX with modifiable plots" is checked
  to_doc <- function(p = NULL, w, h, fun = NULL) {
    if (!isTRUE(modif)) return(invisible(NULL))
    tryCatch({
      if (is.null(fun)) export::graph2ppt(x = p, file = modif_file, append = TRUE,
                                          width = w, height = h)
      else export::graph2ppt(fun = fun, file = modif_file, append = TRUE,
                             width = w, height = h)
    }, error = function(e) cat(file = stderr(), "graph2ppt failed:", e$message, "\n"))
    invisible(NULL)
  }

  ## ---------- 1 : overlap indices ----------
  P  <- !is.na(x)
  P  <- P[rowSums(P) > 0, , drop = FALSE]
  ov <- pairwise_overlap(P)
  pt <- pairs_table(ov, groups)
  gs <- group_summary(P, pt, groups)
  ts <- jaccard_tests(ov$jaccard, groups, nperm = nperm, seed = seed)

  png_hm  <- save_plot(plot_jaccard_heatmap(ov$jaccard, groups,
                         "Pairwise Jaccard index (shared / union of quantified proteins)"),
                       "Jaccard_matrix", 9, 8.5)
  png_ov  <- save_plot(plot_jaccard_heatmap(ov$overlap, groups,
                         "Pairwise overlap coefficient (shared / smaller protein set)"),
                       "Overlap_coefficient_matrix", 9, 8.5)
  png_box <- save_plot(plot_within_group(pt, groups, ts, "Jaccard"),
                       "Jaccard_within_groups", 8, 5)
  png_dsp <- save_plot(plot_dispersion(ts, groups), "Jaccard_dispersion", 7, 5)

  ## ---------- 2 : pvclust ----------
  pv <- run_pvclust(x, nboot = nboot, parallel = parallel, seed = seed)
  et <- pvclust_edge_table(pv, groups)
  png_pv <- save_plot(plot_pvclust_gg(pv, groups, alpha = alpha,
                        title = "Ward's method (ward.D2), correlation-based distance, multiscale bootstrap"),
                      "pvclust_dendrogram", 7.5, 0.28 * ncol(x) + 2)
  # standard pvclust plot (base graphics)
  draw_pv <- function() {
    plot(pv, hang = -1, cex = 0.8, cex.pv = 0.6, main = "Cluster dendrogram with AU/BP values (%)",
         sub = sprintf("Ward.D2, 1 - Pearson correlation (pairwise complete obs.), nboot = %d", nboot),
         xlab = "")
    pvrect(pv, alpha = alpha, pv = "au", border = 2, lty = 2)
  }
  for (dev in c("pdf", "png")) {
    f <- out(paste0("pvclust_standard.", dev))
    if (dev == "pdf") pdf(f, width = 11, height = 7) else png(f, width = 11, height = 7, units = "in", res = 300)
    draw_pv()
    dev.off(); files <- c(files, f)
  }
  png_pvs <- out("pvclust_standard.png")
  to_doc(w = 11, h = 7, fun = draw_pv)

  ## ---------- Excel ----------
  wb <- createWorkbook()
  add <- function(sheet, d, rn = FALSE) {
    addWorksheet(wb, sheet)
    if (rn) d <- data.frame(Sample = rownames(d), d, check.names = FALSE)
    writeDataTable(wb, sheet, d)
  }
  add("README", data.frame(Item = c(
    "Presence", "Jaccard", "Overlap coefficient", "Sorensen",
    "Core_pct_all_samples", "Core_pct_size_matched", "betadisper", "PERMANOVA",
    "pvclust", "AU / BP"), Definition = c(
    "A protein is 'present' in a sample when its intensity is quantified (non-missing) after the app filters.",
    "|A n B| / |A u B| for each pair of samples; independent of the number of samples per group.",
    "|A n B| / min(|A|,|B|); insensitive to differences in the number of proteins identified per sample.",
    "2|A n B| / (|A|+|B|).",
    "% of proteins found in all samples of the group among proteins found in at least one: decreases with group size.",
    sprintf("Same, averaged over all subsets of k = %d samples (k = size of the smallest group) to remove the group-size effect.", min(table(groups))),
    sprintf("Test of homogeneity of multivariate dispersion on the Jaccard distance (1 - J), %d permutations; pairwise tests BH-adjusted.", nperm),
    sprintf("adonis2 on the Jaccard distance, %d permutations: do categories differ in protein composition?", nperm),
    sprintf("hclust Ward.D2 on 1 - Pearson correlation (pairwise complete obs.) of log2 intensities (sample dendrogram of the report); %d bootstrap replicates x 10 scales (r = 0.5..1.4).", nboot),
    "AU: approximately unbiased p-value (multiscale bootstrap); BP: bootstrap probability. Clusters with AU >= 95% are strongly supported.")))
  add("Jaccard_matrix", as.data.frame(ov$jaccard), rn = TRUE)
  add("Overlap_coef_matrix", as.data.frame(ov$overlap), rn = TRUE)
  add("Shared_proteins_matrix", as.data.frame(ov$intersection), rn = TRUE)
  add("Pairs", pt)
  add("Group_summary", gs)
  add("Dispersion_test", ts$dispersion_global)
  add("Dispersion_pairwise", ts$dispersion_pairwise)
  add("Dist_to_centroid", ts$distance_to_centroid)
  add("PERMANOVA", ts$permanova)
  add("pvclust_edges", et)
  add("pvclust_AU95_clusters", {
    pk <- pvpick(pv, alpha = alpha, pv = "au", max.only = FALSE)
    if (length(pk$edges) == 0) data.frame(Edge = NA, Members = "none")
    else et[et$Edge %in% pk$edges, c("Edge", "AU_pct", "BP_pct", "Group_composition", "Members")]
  })
  fx <- out("Jaccard_pvclust_results.xlsx")
  saveWorkbook(wb, fx, overwrite = TRUE)
  files <- c(files, fx)
  saveRDS(list(overlap = ov, pairs = pt, summary = gs, tests = ts[-1], pvclust = pv),
          out("overlap_results.rds"))
  files <- c(files, out("overlap_results.rds"))

  ## ---------- PowerPoint ----------
  if (!is.null(doc)) {
    slide_img <- function(doc, title, png, w, h, txt = NULL) {
      doc <- add_slide(doc, layout = "Title Only", master = "Office Theme")
      doc <- ph_with(doc, fpar(ftext(title, fp_text(font.size = 24))),
                     location = ph_location_type(type = "title"))
      r <- min(9.5 / w, (if (is.null(txt)) 5.6 else 4.9) / h)
      doc <- ph_with(doc, external_img(png, width = w * r, height = h * r),
                     location = ph_location(left = (10 - w * r) / 2, top = 1.5,
                                            width = w * r, height = h * r))
      if (!is.null(txt))
        doc <- ph_with(doc, fpar(ftext(txt, fp_text(font.size = 11))),
                       location = ph_location(left = 0.3, top = 6.5, width = 9.4, height = 0.9))
      doc
    }
    doc <- slide_img(doc, paste0("[", prefix, "] Pairwise Jaccard index"), png_hm, 9, 8.5)
    doc <- slide_img(doc, paste0("[", prefix, "] Pairwise overlap coefficient"), png_ov, 9, 8.5)
    doc <- slide_img(doc, paste0("[", prefix, "] Within-group Jaccard index"), png_box, 8, 5,
                     txt = paste(sprintf("%s: median J = %.2f [%.2f-%.2f]", gs$Group,
                                         gs$Jaccard_within_median, gs$Jaccard_within_Q1,
                                         gs$Jaccard_within_Q3), collapse = "; "))
    doc <- slide_img(doc, paste0("[", prefix, "] Within-group heterogeneity"), png_dsp, 7, 5)
    doc <- slide_img(doc, paste0("[", prefix, "] Clustering with AU/BP p-values (pvclust)"),
                     png_pv, 7.5, 0.28 * ncol(x) + 2)
    doc <- slide_img(doc, paste0("[", prefix, "] pvclust dendrogram (standard plot)"),
                     png_pvs, 11, 7)
  }

  list(doc = doc, files = files,
       results = list(overlap = ov, pairs = pt, summary = gs, tests = ts, pvclust = pv, edges = et))
}
