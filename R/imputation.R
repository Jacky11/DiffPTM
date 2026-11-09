#' Impute missing values
#'
#' Wraps `imp4p` (Giai Gianetto et al.) as in the original app. `"none"` leaves `NA`s in
#' place; the limma method tolerates them (row-wise fitting), which is the recommended
#' setting for the DPU analysis where imputed PTM and protein values would otherwise
#' propagate into the corrected fold-change.
#'
#' `raw_detected` is never modified, so detection-based analyses
#' ([detection_discordance()], [qc_jaccard()]) always use the real observation pattern.
#'
#' @param x a `ptm_experiment`.
#' @param assay `"protein"`, `"ptm"` or both.
#' @param method `"none"`, `"mle"` (imp4p::impute.mle), `"slsa"` (imp4p::impute.slsa),
#'   `"mi"` (imp4p::impute.mi, MCAR/MNAR mixture), `"min_det"` (per-sample minimum, quick
#'   placeholder for development).
#' @export
impute <- function(x, assay = c("protein", "ptm"), method = c("none", "mle", "slsa", "mi", "min_det")) {
  method <- match.arg(method)
  if (method == "none") return(add_log(x, "No imputation; missing values handled by the statistical model."))
  for (a in assay) {
    if (is.null(x[[a]])) next
    M <- x[[a]]$intensity
    g <- groups_of(x, a)
    if (method %in% c("mle", "slsa", "mi") && !requireNamespace("imp4p", quietly = TRUE))
      stop("Package imp4p is required for method ", method)
    M2 <- switch(method,
      mle     = imp4p::impute.mle(M, conditions = g),
      slsa    = imp4p::impute.slsa(M, conditions = g),
      mi      = imp4p::impute.mi(M, conditions = g),
      min_det = { mins <- apply(M, 2, min, na.rm = TRUE); for (j in seq_len(ncol(M))) M[is.na(M[, j]), j] <- mins[j]; M }
    )
    M2 <- as.matrix(M2); dimnames(M2) <- dimnames(M)
    x[[a]]$intensity <- M2
    x <- add_log(x, sprintf("%s: imputation of %d missing values with imp4p::%s.", a, sum(is.na(M)),
                            switch(method, mle = "impute.mle", slsa = "impute.slsa", mi = "impute.mi", min_det = "min_det")))
  }
  x
}
