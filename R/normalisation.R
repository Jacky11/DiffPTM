#' Normalise an assay
#'
#' Each assay (sites, proteins) is centred separately, with one factor per sample applied to
#' every feature of the assay.
#'
#' * `"global_median"` (default): median centering across all samples, as MSstats'
#'   `equalizeMedians` and MSstatsTMT's global normalisation. A shift of a whole assay between
#'   conditions (a difference of loading or enrichment efficiency, or a genuine global change)
#'   is removed.
#' * `"median"`: per-condition median centering (the default of DiffPTM <= 2.0.0) — within each
#'   group, sample medians are computed on features without missing values in that group and
#'   each sample is shifted to the mean of the group medians, so that differences between groups
#'   are preserved, including a difference of mean loading between conditions. For the usage
#'   test, such a difference in the PTM assay is attributed to every site: in simulations with
#'   random per-sample loading offsets (SD 0.3 log2) of the PTM assay, it raised the mean FDP of
#'   every strategy to about 0.2 with three samples per condition, whereas global centering kept
#'   it at the nominal level. Use it only when a global change of the enriched fraction is
#'   expected to be biological and loading is controlled. Groups with no complete feature are
#'   centred on all available values and reported in the log.
#' * `"combat"`: `sva::ComBat` with the `factor` column (batch) of the annotation.
#' * `"none"`.
#'
#' @param x a `ptm_experiment`.
#' @param assay `"protein"`, `"ptm"` or both.
#' @param method see above.
#' @export
normalise <- function(x, assay = c("protein", "ptm"), method = c("global_median", "median", "combat", "none")) {
  method <- match.arg(method)
  for (a in assay) {
    if (is.null(x[[a]]) || method == "none") next
    M <- x[[a]]$intensity
    g <- groups_of(x, a)
    if (method == "median") {
      out <- median_norm_by_group(M, g)
      M <- out$intensity
      note <- if (length(out$fallback)) sprintf(" (groups without complete features, centred on all values: %s)", paste(out$fallback, collapse = ", ")) else ""
      x <- add_log(x, sprintf("%s: per-condition median centering%s.", a, note))
    } else if (method == "global_median") {
      med <- apply(M, 2, stats::median, na.rm = TRUE)
      M <- sweep(M, 2, med - mean(med))
      x <- add_log(x, sprintf("%s: global median centering.", a))
    } else if (method == "combat") {
      if (!requireNamespace("sva", quietly = TRUE)) stop("Package sva is required for ComBat")
      ann <- x$annotation[[a]]
      if (all(is.na(ann$factor))) stop("ComBat needs the `factor` column of the annotation as batch")
      mod <- stats::model.matrix(~ factor(ann$condition))
      M <- sva::ComBat(dat = M, batch = factor(ann$factor), mod = mod)
      x <- add_log(x, sprintf("%s: ComBat batch correction using `%s` as batch.", a, "factor"))
    }
    x[[a]]$intensity <- M
  }
  x
}

median_norm_by_group <- function(M, g) {
  out <- M; fallback <- character(0)
  for (lv in levels(g)) {
    cols <- which(g == lv)
    sub <- M[, cols, drop = FALSE]
    complete <- stats::complete.cases(sub)
    if (!any(complete)) { complete <- rep(TRUE, nrow(sub)); fallback <- c(fallback, lv) }
    med <- apply(sub[complete, , drop = FALSE], 2, stats::median, na.rm = TRUE)
    out[, cols] <- sweep(sub, 2, med - mean(med))
  }
  list(intensity = out, fallback = fallback)
}
