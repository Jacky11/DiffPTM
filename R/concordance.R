#' Compare the results of several methods
#'
#' @param res tidy table from [run_methods()] with >= 2 methods.
#' @param fdr,min_log2fc significance thresholds.
#' @return list with `wide` (one row per feature x analysis x contrast, one `log2FC_<method>` /
#'   `sig_<method>` column per method), `overlap` (counts of significant features per method
#'   combination) and `correlation` (Pearson correlation of log2FC between methods).
#' @export
compare_methods <- function(res, fdr = 0.01, min_log2fc = 1) {
  res$sig <- !is.na(res$adj_pvalue) & res$adj_pvalue < fdr & abs(res$log2FC) >= min_log2fc
  methods <- unique(res$method)
  key <- unique(res[, c("feature_id", "analysis", "contrast")])
  wide <- key
  for (e in methods) {
    sub <- res[res$method == e, ]
    idx <- match(paste(key$feature_id, key$analysis, key$contrast), paste(sub$feature_id, sub$analysis, sub$contrast))
    wide[[paste0("log2FC_", e)]] <- sub$log2FC[idx]
    wide[[paste0("adj_pvalue_", e)]] <- sub$adj_pvalue[idx]
    wide[[paste0("sig_", e)]] <- sub$sig[idx] %in% TRUE
  }
  sig_cols <- paste0("sig_", methods)
  wide$n_methods_significant <- rowSums(wide[, sig_cols, drop = FALSE])
  pattern <- apply(wide[, sig_cols, drop = FALSE], 1, function(v) paste(methods[v], collapse = "&"))
  pattern[pattern == ""] <- "none"
  overlap <- as.data.frame(table(analysis = wide$analysis, contrast = wide$contrast, methods = pattern))
  overlap <- overlap[overlap$Freq > 0, ]
  correlation <- NULL
  if (length(methods) >= 2) {
    cmb <- utils::combn(methods, 2)
    correlation <- do.call(rbind, lapply(split(wide, paste(wide$analysis, wide$contrast)), function(w) {
      do.call(rbind, lapply(seq_len(ncol(cmb)), function(k) {
        a <- w[[paste0("log2FC_", cmb[1, k])]]; b <- w[[paste0("log2FC_", cmb[2, k])]]
        data.frame(analysis = w$analysis[1], contrast = w$contrast[1], method_a = cmb[1, k], method_b = cmb[2, k],
                   n = sum(!is.na(a) & !is.na(b)), pearson = suppressWarnings(stats::cor(a, b, use = "complete.obs")))
      }))
    }))
    rownames(correlation) <- NULL
  }
  list(wide = wide, overlap = overlap, correlation = correlation)
}

#' Consensus calls: features significant in at least `min_methods` methods
#' @inheritParams compare_methods
#' @param min_methods minimal number of methods calling the feature.
#' @export
consensus_calls <- function(res, fdr = 0.01, min_log2fc = 1, min_methods = 2) {
  w <- compare_methods(res, fdr, min_log2fc)$wide
  w[w$n_methods_significant >= min_methods, ]
}

#' UpSet plot of the features called by each combination of methods
#'
#' Bars give the number of features significant with exactly the methods marked by the dots below
#' (an UpSet plot, which reads like a Venn diagram for any number of methods); the number of
#' features called by each method is given next to its name.
#' @inheritParams compare_methods
#' @param analysis,contrast the analysis (`"ptm_DPU"`, `"ptm_DPA"`, `"protein_DPA"`) and contrast shown.
#' @param title plot title (default: analysis and contrast).
#' @param subtitle plot subtitle (`NULL` for none).
#' @return a ggplot object, or `NULL` when fewer than two methods were run.
#' @export
plot_method_overlap <- function(res, fdr = 0.05, min_log2fc = 0, analysis = "ptm_DPU", contrast = NULL, title = NULL,
                                subtitle = "Features significant with exactly the marked methods; (n) = calls per method") {
  methods <- unique(res$method)
  if (length(methods) < 2) return(NULL)
  r <- res[res$analysis == analysis, ]
  if (is.null(contrast)) contrast <- unique(r$contrast)[1]
  r <- r[r$contrast == contrast, ]
  w <- compare_methods(r, fdr, min_log2fc)$wide
  S <- as.matrix(w[, paste0("sig_", methods), drop = FALSE]); colnames(S) <- methods
  if (is.null(title)) title <- sprintf("Agreement between methods: %s, %s", analysis, contrast)
  if (!any(S)) return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = "No significant feature with any method") +
                        ggplot2::theme_void() + ggplot2::labs(title = title))
  plot_upset(S, title, subtitle)
}
