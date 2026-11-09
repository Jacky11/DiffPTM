#' Remove flagged features (reverse hits, contaminants, "only identified by site")
#'
#' Flag columns are created by the importers from the software's flag columns
#' (logical columns named after `make.names()` of the original column).
#'
#' @param x a `ptm_experiment`.
#' @param assay `"protein"`, `"ptm"` or both.
#' @param flags names of logical columns in `meta`; `NULL` = every logical column whose
#'   name matches `Reverse|contaminant|Only.identified`.
#' @export
filter_flagged <- function(x, assay = c("protein", "ptm"), flags = NULL) {
  for (a in assay) {
    if (is.null(x[[a]])) next
    meta <- x[[a]]$meta
    fl <- flags %||% grep("Reverse|contaminant|Only\\.identified", names(meta), value = TRUE, ignore.case = TRUE)
    fl <- fl[fl %in% names(meta)]
    if (!length(fl)) next
    drop <- Reduce(`|`, lapply(fl, function(f) isTRUE_vec(meta[[f]])))
    x[[a]] <- subset_assay(x[[a]], !drop)
    x <- add_log(x, sprintf("%s: removed %d features flagged in %s.", a, sum(drop), paste(fl, collapse = ", ")))
  }
  x
}

isTRUE_vec <- function(v) !is.na(v) & as.logical(v)

#' Filter PTM features on localization probability
#'
#' @param x a `ptm_experiment`.
#' @param keep_na keep features whose probability is `NA`.
#' @param min_prob minimum probability (features with `NA` probability are kept when
#'   `keep_na = TRUE`, e.g. Proteome Discoverer tables without a probability column).
#' @export
filter_localization <- function(x, min_prob = 0.75, keep_na = TRUE) {
  if (is.null(x$ptm)) return(x)
  p <- x$ptm$meta$localization_prob
  keep <- if (keep_na) is.na(p) | p >= min_prob else !is.na(p) & p >= min_prob
  x$ptm <- subset_assay(x$ptm, keep)
  add_log(x, sprintf("ptm: kept %d/%d features with localization probability >= %s.",
                     sum(keep), length(keep), min_prob))
}

#' Filter on the number / percentage of observed values
#'
#' A feature is kept when at least one group satisfies both `min_n` observed values and
#' `min_pct` percent observed values (same rule as the original app, vectorised).
#'
#' @inheritParams filter_flagged
#' @param min_n minimal number of observed values in at least one group.
#' @param min_pct minimal percentage (0-100) of observed values in at least one group.
#' @param in_all_groups if `TRUE` the rule must hold in every group instead of at least one.
#' @export
filter_min_observed <- function(x, min_n = 2, min_pct = 0, assay = c("protein", "ptm"), in_all_groups = FALSE) {
  if (is.character(min_n)) stop("`min_n` must be a number; pass the assay as `assay = `")
  for (a in assay) {
    if (is.null(x[[a]])) next
    g <- groups_of(x, a)
    obs <- !is.na(x[[a]]$intensity)
    pass <- sapply(levels(g), function(lv) {
      cols <- which(g == lv)
      n <- rowSums(obs[, cols, drop = FALSE])
      n >= min_n & (100 * n / length(cols)) >= min_pct
    })
    if (is.null(dim(pass))) pass <- matrix(pass, ncol = nlevels(g))
    keep <- if (in_all_groups) rowSums(pass) == ncol(pass) else rowSums(pass) > 0
    x[[a]] <- subset_assay(x[[a]], keep)
    x <- add_log(x, sprintf("%s: kept %d/%d features with >= %d observed values and >= %s%% observed values in %s group.",
                            a, sum(keep), length(keep), min_n, min_pct, if (in_all_groups) "every" else "at least one"))
  }
  x
}

#' Subset an assay by a logical / index vector on features
#' @keywords internal
subset_assay <- function(assay, keep) {
  assay$intensity    <- assay$intensity[keep, , drop = FALSE]
  assay$raw_detected <- assay$raw_detected[keep, , drop = FALSE]
  assay$meta         <- assay$meta[keep, , drop = FALSE]
  for (m in c("protein_intensity", "n_features", "protein_n_features"))
    if (!is.null(assay[[m]])) assay[[m]] <- assay[[m]][keep, , drop = FALSE]
  rownames(assay$meta) <- NULL
  assay
}

#' Groups with no complete feature (cannot be normalised by median of complete rows)
#' @keywords internal
groups_without_complete_rows <- function(mat, g) {
  vapply(levels(g), function(lv) {
    sub <- mat[, g == lv, drop = FALSE]
    !any(stats::complete.cases(sub))
  }, logical(1))
}
