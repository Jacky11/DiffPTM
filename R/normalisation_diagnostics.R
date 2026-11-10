#' Help to choose between global and per-condition centering
#'
#' Run on filtered, *not yet normalised* data. For each assay and each pair of groups that the
#' contrasts compare (by default every pair of conditions, within each level of the second factor
#' if there is one, see [make_contrasts()]), it reports:
#'
#' * `shift`: the difference between the groups of the mean sample medians (features observed
#'   in every sample of the two groups), with a t-test on the sample medians (`shift_p`): is
#'   there a shift between conditions beyond the variation of the samples within conditions?
#' * `median_lfc` and `bulk_lfc`: the median and the mode (most common value, from a kernel
#'   density) of the per-feature log2 fold changes. A loading or enrichment difference moves
#'   every feature, hence the bulk (`bulk_lfc` close to `median_lfc`); a change of a subset of
#'   features moves the median but leaves the bulk in place.
#'
#' The recommendation, made on the PTM assay, is:
#'
#' 1. no shift (not significant, or smaller than `min_shift` log2): **global** centering (both
#'    choices give nearly the same results);
#' 2. shift, and the bulk did not move (`|bulk_lfc| < bulk_ratio * |median_lfc|`): the
#'    shift comes from a subset of sites, which global centering would spread to every site:
#'    **per-condition** centering;
#' 3. shift of the whole distribution: a loading or enrichment difference, or a
#'    genuine global change of the modification; the data cannot tell them apart, so
#'    **global** centering unless loading and enrichment were controlled independently (equal
#'    input checked on the total proteome, spiked standards, reference sites). The message says
#'    whether the bulk of the proteome of the same samples moved in the same direction by at least
#'    half as much, which points to loading.
#'
#' The p-values are adjusted over the pairs of each assay (Holm). This is a diagnostic: the choice
#' remains the user's, and it is reported in the Methods paragraph. The rule was evaluated in
#' simulations with technical shifts, loading differences and biological changes of a subset or of
#' all sites (see the DiffPTM article, Supplementary Note 3).
#' @param x a `ptm_experiment`, filtered but not normalised.
#' @param contrasts contrasts as returned by [make_contrasts()] (`group1`, `group2`); default: every
#'   pair of conditions (within each level of the second factor).
#' @param alpha level of the test of the shift.
#' @param bulk_ratio the bulk is considered not to have moved when the mode of the fold changes
#'   is smaller than this fraction of their median (in absolute value).
#' @param min_shift minimal shift (log2) between conditions considered (smaller shifts barely
#'   change the results).
#' @param min_features minimal number of features observed in every sample of two groups.
#' @return a list with `table` (one row per assay and compared pair of groups), `recommendation` (`"global_median"`
#'   or `"median"`), `case` (1, 2 or 3, see above) and `message` (plain-language explanation).
#' @examples
#' sim <- simulate_ptm_features(n_rep = 3, n_proteins = 300, n_sites = 250, seed = 1)
#' x <- sim_to_experiment(sim)
#' normalisation_diagnostics(x)$message
#' @export
normalisation_diagnostics <- function(x, contrasts = NULL, alpha = 0.05, bulk_ratio = 0.5, min_shift = 0.1, min_features = 50) {
  if (any(grepl("median centering|ComBat", x$log))) warning("x seems already normalised: run normalisation_diagnostics() before normalise()")
  mode_of <- function(v) { d <- stats::density(v, n = 4096); d$x[which.max(d$y)] }
  rows <- list()
  for (a in intersect(c("ptm", "protein"), names(x))) {
    if (is.null(x[[a]])) next
    g <- groups_of(x, a); M <- x[[a]]$intensity; lv <- levels(g)
    if (length(lv) < 2) next
    ctr <- contrasts %||% tryCatch(make_contrasts(x, "pairwise"), error = function(e) NULL)
    if (is.null(ctr)) next
    ctr <- ctr[ctr$group1 %in% lv & ctr$group2 %in% lv, , drop = FALSE]                 # vs_all contrasts are skipped
    for (i in seq_len(nrow(ctr))) {
      k <- ctr$group1[i]; r0 <- ctr$group2[i]
      j1 <- which(g == k); j0 <- which(g == r0)
      cc <- stats::complete.cases(M[, c(j0, j1), drop = FALSE])
      if (sum(cc) < min_features || length(j0) < 2 || length(j1) < 2) {
        rows[[length(rows) + 1]] <- data.frame(assay = a, group = k, reference = r0, n_features = sum(cc), shift = NA_real_,
                                               shift_p = NA_real_, median_lfc = NA_real_, bulk_lfc = NA_real_)
        next
      }
      med <- apply(M[cc, c(j0, j1), drop = FALSE], 2, stats::median)
      m0 <- med[seq_along(j0)]; m1 <- med[-seq_along(j0)]
      p <- tryCatch(stats::t.test(m1, m0, var.equal = TRUE)$p.value, error = function(e) NA_real_)
      lfc <- rowMeans(M[cc, j1, drop = FALSE]) - rowMeans(M[cc, j0, drop = FALSE])
      rows[[length(rows) + 1]] <- data.frame(assay = a, group = k, reference = r0, n_features = sum(cc), shift = mean(m1) - mean(m0),
                                             shift_p = p, median_lfc = stats::median(lfc), bulk_lfc = mode_of(lfc))
    }
  }
  tab <- do.call(rbind, rows)
  if (is.null(tab)) return(list(table = NULL, recommendation = "global_median", case = NA_integer_,
                                message = "Fewer than two groups: global median centering."))
  tab$shift_p_adj <- stats::ave(tab$shift_p, tab$assay, FUN = function(p) stats::p.adjust(p, "holm"))
  tab$significant <- !is.na(tab$shift_p_adj) & tab$shift_p_adj < alpha & abs(tab$shift) >= min_shift
  tab$bulk_moved <- !is.na(tab$bulk_lfc) & abs(tab$bulk_lfc) >= bulk_ratio * abs(tab$median_lfc)
  pt <- tab[tab$assay == "ptm", , drop = FALSE]
  if (!nrow(pt)) pt <- tab
  f2 <- function(v) sprintf("%+.2f", v)
  sig <- pt[pt$significant, , drop = FALSE]
  if (!nrow(sig)) {
    case <- 1L; rec <- "global_median"
    msg <- sprintf("No shift of the PTM assay between conditions (largest: %s log2, p = %.2g; a shift counts when significant and at least %.2g log2): both centerings give nearly the same results; global median centering (default).",
                   f2(pt$shift[which.max(abs(pt$shift))]), pt$shift_p_adj[which.max(abs(pt$shift))], min_shift)
  } else if (!any(sig$bulk_moved)) {
    case <- 2L; rec <- "median"
    msg <- sprintf("The PTM assay is shifted between conditions (%s), but the bulk of the sites did not move (mode of the fold changes %s): the shift comes from a subset of sites, which global centering would spread to every site. Per-condition median centering is recommended.",
                   paste(sprintf("%s vs %s: median fold change %s, p = %.2g", sig$group, sig$reference, f2(sig$median_lfc), sig$shift_p_adj), collapse = "; "),
                   paste(f2(sig$bulk_lfc), collapse = ", "))
  } else {
    case <- 3L; rec <- "global_median"
    b <- sig[sig$bulk_moved, , drop = FALSE]
    pr <- tab[tab$assay == "protein", , drop = FALSE]
    k <- match(paste(b$group, b$reference), paste(pr$group, pr$reference))
    shared <- !is.na(k) & !is.na(pr$bulk_lfc[k]) & sign(pr$bulk_lfc[k]) == sign(b$bulk_lfc) & abs(pr$bulk_lfc[k]) >= 0.5 * abs(b$bulk_lfc)   # loading moves both fractions alike
    pairs <- function(v) paste(sprintf("%s vs %s", b$group[v], b$reference[v]), collapse = ", ")
    msg <- paste0(sprintf("The whole PTM distribution is shifted between conditions (%s). ",
                          paste(sprintf("%s vs %s: %s log2, adjusted p = %.2g", b$group, b$reference, f2(b$bulk_lfc), b$shift_p_adj), collapse = "; ")),
                  if (all(shared)) "The proteome of the same samples is shifted in the same direction, which points to a loading difference. "
                  else if (any(shared)) sprintf("The proteome of the same samples is shifted in the same direction for %s (pointing to a loading difference), not for %s. ", pairs(shared), pairs(!shared))
                  else "The proteome is not shifted in the same way: a difference of enrichment efficiency or a genuine global change of the modification. ",
                  "The data cannot tell a technical shift from a global biological change: global median centering is recommended, unless loading and enrichment were controlled independently (spiked standards, reference sites), in which case per-condition centering keeps a global change.")
  }
  list(table = tab, recommendation = rec, case = case, message = msg)
}
