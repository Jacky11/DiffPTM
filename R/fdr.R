#' Adjust p-values with a pi0 estimate (cp4p) or plain Benjamini-Hochberg
#'
#' @param p numeric vector of p-values (NA allowed).
#' @param pi0_method one of the `cp4p::estim.pi0` methods (`"abh"`, `"st.spline"`,
#'   `"st.boot"`, `"langaas"`, `"histo"`, `"pounds"`, `"jiang"`, `"slim"`) or `"bh"`.
#' @return list(adjusted, pi0, method).
#' @export
adjust_pvalues <- function(p, pi0_method = "pounds") {
  ok <- !is.na(p)
  adj <- rep(NA_real_, length(p)); pi0 <- 1
  if (pi0_method == "bh" || !requireNamespace("cp4p", quietly = TRUE)) {
    adj[ok] <- stats::p.adjust(p[ok], method = "BH")
    return(list(adjusted = adj, pi0 = 1, method = "bh"))
  }
  res <- cp4p::adjust.p(p[ok], pi0.method = pi0_method)
  adj[ok] <- res$adjp$adjusted.p
  list(adjusted = adj, pi0 = res$pi0, method = pi0_method)
}
