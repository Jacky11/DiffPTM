#' Build the list of contrasts between groups
#'
#' @param x a `ptm_experiment` (groups are taken from the PTM annotation when present,
#'   otherwise the protein annotation; both assays must use the same group labels).
#' @param type `"pairwise"` (every pair of conditions), `"vs_reference"` (each condition
#'   against `reference`, the usual design with a control group) or `"vs_all"` (each
#'   condition against the average of the others). With a second factor, pairs are built
#'   within each factor level (as the original app) unless `within_factor = FALSE`.
#' @param reference the control condition for `"vs_reference"` (default: the first group).
#' @param within_factor build contrasts within each level of the second factor.
#' @return data.frame(name, group1, group2) — for `"vs_all"`, `group2` is a `+`-separated
#'   list of groups whose average is used.
#' @export
make_contrasts <- function(x, type = c("pairwise", "vs_reference", "vs_all"), reference = NULL, within_factor = TRUE) {
  type <- match.arg(type)
  a <- if (!is.null(x$annotation$ptm)) "ptm" else "protein"
  ann <- x$annotation[[a]]
  g <- levels(groups_of(x, a))
  has_factor <- !all(is.na(ann$factor))
  blocks <- if (has_factor && within_factor) split(g, sub("^.*\\.", "", g)) else list(all = g)
  rows <- list()
  for (b in blocks) {
    if (type == "pairwise") {
      if (length(b) < 2) next
      cmb <- utils::combn(b, 2)
      for (k in seq_len(ncol(cmb)))
        rows[[length(rows) + 1]] <- data.frame(name = paste0(cmb[1, k], "_vs_", cmb[2, k]),
                                               group1 = cmb[1, k], group2 = cmb[2, k], stringsAsFactors = FALSE)
    } else if (type == "vs_reference") {
      ref <- if (is.null(reference) || !nzchar(reference)) b[1] else {
        hit <- b[b == reference | sub("\\..*$", "", b) == reference]
        if (!length(hit)) stop("Reference condition not found: ", reference); hit[1] }
      for (gi in setdiff(b, ref))
        rows[[length(rows) + 1]] <- data.frame(name = paste0(gi, "_vs_", ref), group1 = gi, group2 = ref, stringsAsFactors = FALSE)
    } else {
      for (gi in b)
        rows[[length(rows) + 1]] <- data.frame(name = paste0(gi, "_vs_others"), group1 = gi,
                                               group2 = paste(setdiff(b, gi), collapse = "+"), stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) stop("No contrast could be built (need >= 2 groups)")
  do.call(rbind, rows)
}

#' Design matrix without intercept (one column per group), optional subject blocking
#'
#' @param group factor of group labels.
#' @param subject optional subject factor (fixed-effect blocking).
#' @param mixture optional TMT plex factor (fixed batch effect); ignored when `subject` is
#'   given, since samples are nested in plexes.
#' @export
make_design <- function(group, subject = NULL, mixture = NULL) {
  group <- factor(group)
  lv <- levels(group)
  use_subj <- !is.null(subject) && !all(is.na(subject))
  use_mix <- !is.null(mixture) && !all(is.na(mixture)) && length(unique(mixture)) > 1
  d <- data.frame(group = group)
  f <- "~ 0 + group"
  if (use_subj) { d$subject <- factor(subject); f <- paste(f, "+ subject") }
  if (use_mix)  { d$mixture <- factor(mixture); f <- paste(f, "+ mixture") }
  design <- stats::model.matrix(stats::as.formula(f), data = d)
  # drop columns that are linear combinations of earlier ones (e.g. a mixture effect already carried by
  # the sample terms); group columns come first and are always kept
  q <- qr(design)
  if (q$rank < ncol(design)) design <- design[, sort(q$pivot[seq_len(q$rank)]), drop = FALSE]
  colnames(design) <- sub("^group", "", colnames(design))
  colnames(design)[seq_along(lv)] <- make.names(lv)
  colnames(design)[-seq_along(lv)] <- make.names(colnames(design)[-seq_along(lv)], unique = TRUE)
  attr(design, "group_levels") <- stats::setNames(make.names(lv), lv)
  design
}

#' Methods currently usable (installed dependencies)
#' @export
available_methods <- function() {
  e <- c(limma = TRUE,
         msstatsptm = requireNamespace("MSstatsPTM", quietly = TRUE),
         msqrob2ptm = requireNamespace("msqrob2", quietly = TRUE) && requireNamespace("QFeatures", quietly = TRUE))
  names(e)[e]
}

#' Run one or several statistical methods and return one tidy table
#'
#' @param x a `ptm_experiment` after filtering/normalisation (and matching for DPU).
#' @param contrasts from [make_contrasts()].
#' @param methods subset of `c("limma", "msstatsptm", "msqrob2ptm")`.
#' @param analyses subset of `c("protein_DPA", "ptm_DPA", "ptm_DPU")`.
#' @param pi0_method see [adjust_pvalues()] (used by the limma method; the other methods
#'   return their own adjusted p-values, re-adjusted only if `readjust = TRUE`).
#' @param readjust re-adjust the p-values of the other methods with `pi0_method`.
#' @param ... passed to the methods.
#' @return data.frame with columns `feature_id, analysis, contrast, log2FC, se, statistic,
#'   pvalue, adj_pvalue, method`.
#' @export
run_methods <- function(x, contrasts, methods = "limma",
                        analyses = c("protein_DPA", "ptm_DPA", "ptm_DPU"),
                        pi0_method = "pounds", readjust = FALSE, ...) {
  methods <- match.arg(methods, c("limma", "msstatsptm", "msqrob2ptm"), several.ok = TRUE)
  if (is.null(contrasts) || !nrow(contrasts)) stop("No contrast to test")
  out <- lapply(methods, function(e) {
    cli::cli_alert_info("Running method {e}")
    res <- switch(e,
      limma      = method_limma(x, contrasts, analyses = analyses, pi0_method = pi0_method, ...),
      msstatsptm = method_msstatsptm(x, contrasts, analyses = analyses, ...),
      msqrob2ptm = method_msqrob2ptm(x, contrasts, analyses = analyses, ...))
    if (readjust && e != "limma") {
      for (k in split(seq_len(nrow(res)), paste(res$analysis, res$contrast)))
        res$adj_pvalue[k] <- adjust_pvalues(res$pvalue[k], pi0_method)$adjusted
    }
    res
  })
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  add_confidence_interval(res)
}

#' Empty result table with the common schema
#' @keywords internal
empty_result <- function() {
  data.frame(feature_id = character(0), analysis = character(0), contrast = character(0),
             log2FC = numeric(0), se = numeric(0), df = numeric(0), statistic = numeric(0), pvalue = numeric(0),
             adj_pvalue = numeric(0), ci_lo = numeric(0), ci_hi = numeric(0), method = character(0), stringsAsFactors = FALSE)
}

#' Add a 95 % confidence interval to a result table
#'
#' `ci_lo` / `ci_hi` = `log2FC` ± t(0.975, `df`) × `se`, using the degrees of freedom each method
#' reports (moderated degrees of freedom for limma and msqrob2PTM, Satterthwaite for MSstatsPTM);
#' `NA` when `df` or `se` is missing. The interval is that of the method's own t-statistic, so it is
#' consistent with its p-value.
#' @param res a result table.
#' @param level confidence level.
#' @return `res` with `ci_lo` and `ci_hi` columns (placed after `adj_pvalue`).
#' @export
add_confidence_interval <- function(res, level = 0.95) {
  df <- if ("df" %in% names(res)) res$df else NA_real_
  q <- stats::qt(1 - (1 - level) / 2, df)
  res$ci_lo <- res$log2FC - q * res$se; res$ci_hi <- res$log2FC + q * res$se
  first <- c("feature_id", "analysis", "contrast", "log2FC", "se", "df", "statistic", "pvalue", "adj_pvalue", "ci_lo", "ci_hi", "method")
  res[, c(intersect(first, names(res)), setdiff(names(res), first)), drop = FALSE]
}
