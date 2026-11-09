#' MSstatsPTM method ("model then correct"), harmonised inputs
#'
#' Adapter around `MSstatsPTM::dataSummarizationPTM()` (label-free) or
#' `MSstatsPTM::dataSummarizationPTM_TMT()` (TMT) and `MSstatsPTM::groupComparisonPTM()`,
#' run on the **same processed values as the other methods**: one feature per site and per
#' protein (DiffPTM's log2 intensities, passed as `2^log2`), with MSstats' normalisation and
#' imputation turned off. This isolates the correction strategy. To run MSstatsPTM as its
#' users do (feature-level input, default normalisation and AFT imputation), use
#' [method_msstatsptm_native()].
#'
#' PTM features are named `<protein_key>_<site_key>` with opaque keys (`P12_S34`): MSstatsPTM
#' finds the global protein of a site from that name, and real identifiers may contain `_`,
#' `;` or `.` (protein groups), which breaks the match. `ptm_DPU` is MSstatsPTM's
#' `ADJUSTED.Model` restricted to `Adjusted == TRUE` rows, with BH re-applied to them;
#' `ptm_DPA` is `PTM.Model`; `protein_DPA` is `PROTEIN.Model`.
#'
#' TMT: the annotation's `mixture` column gives the plex (`Mixture`); the sample is the
#' `Channel`; MSstatsTMT reference-channel normalisation is off (no reference channel).
#'
#' @inheritParams run_methods
#' @param ... ignored.
#' @param data_type `"LabelFree"` or `"TMT"`.
#' @param moderated use MSstats' empirical-Bayes moderated t-test (default `FALSE`, as in MSstatsPTM).
#' @export
method_msstatsptm <- function(x, contrasts, analyses = c("protein_DPA", "ptm_DPA", "ptm_DPU"),
                              data_type = c("LabelFree", "TMT"), moderated = FALSE, ...) {
  data_type <- match.arg(data_type)
  if (!requireNamespace("MSstatsPTM", quietly = TRUE)) stop("Package MSstatsPTM is not installed")
  if (is.null(x$ptm$protein_intensity)) stop("Run match_ptm_to_protein() first")

  parent_ids <- x$ptm$meta$matched_protein                          # "" when unmatched
  prot_all_ids <- x$protein$meta$feature_id
  pseudo_ids <- setdiff(unique(parent_ids[nzchar(parent_ids)]), prot_all_ids)   # "a;b": mean of parents
  prot_map <- data.frame(feature_id = c(prot_all_ids, pseudo_ids), stringsAsFactors = FALSE)
  prot_map$key <- paste0("P", seq_len(nrow(prot_map)))
  prot_key <- prot_map$key[match(parent_ids, prot_map$feature_id)]
  prot_key[is.na(prot_key)] <- "P0"
  ptm_name <- paste0(prot_key, "_S", seq_len(nrow(x$ptm$meta)))
  site_map <- data.frame(key = ptm_name, feature_id = x$ptm$meta$feature_id, stringsAsFactors = FALSE)

  # PROTEIN dataset = the whole proteome plus one pseudo-protein per multi-parent key
  P <- rbind(x$protein$intensity, x$ptm$protein_intensity[match(pseudo_ids, parent_ids), , drop = FALSE])
  to_long <- if (data_type == "TMT") msstats_long_tmt else msstats_long_lf
  input <- list(PTM = to_long(x$ptm$intensity, ptm_name, x$annotation$ptm),
                PROTEIN = to_long(P, prot_map$key, x$annotation$protein))
  settings <- if (data_type == "TMT")
    list(method = "msstats", global_norm = FALSE, global_norm.PTM = FALSE, reference_norm = FALSE,
         reference_norm.PTM = FALSE, MBimpute = FALSE, MBimpute.PTM = FALSE)
  else list(normalization = FALSE, normalization.PTM = FALSE, MBimpute = FALSE, MBimpute.PTM = FALSE)
  gc <- msstats_fit(input, contrasts, data_type, settings, moderated)
  res <- msstats_tidy(gc, analyses, site_map,
                      prot_map = data.frame(key = prot_map$key, feature_id = prot_map$feature_id),
                      keep_proteins = prot_all_ids, method = "msstatsptm")
  res
}

#' MSstatsPTM method as its users run it ("native")
#'
#' Runs MSstatsPTM from feature-level input (peptide x charge, or PSM for TMT) with its own
#' processing: median-equalisation normalisation, Tukey median polish summarisation and,
#' by default, accelerated-failure-time imputation of censored values (`MBimpute = TRUE`).
#' This is the reference arm for the comparison of complete pipelines, alongside the
#' harmonised [method_msstatsptm()].
#'
#' @param input list(PTM = , PROTEIN = ) in MSstats format, e.g. from
#'   `MSstatsPTM::MaxQtoMSstatsPTMFormat()` / `MSstats::MaxQtoMSstatsFormat()` (see the
#'   MSstatsPTM tutorial) or [sim_to_msstats_input()]. PTM `ProteinName` must be
#'   `<protein>_<site>` with `<protein>` equal to the PROTEIN `ProteinName`.
#' @param contrasts data.frame(name, group1, group2) whose groups are values of the
#'   `Condition` column (e.g. from [make_contrasts()]).
#' @param data_type `"LabelFree"` or `"TMT"`.
#' @param MBimpute MSstats' AFT imputation (default `TRUE`, MSstats' default).
#' @param normalization `"equalizeMedians"` (default) or `FALSE`; for TMT, global
#'   normalisation on/off.
#' @param moderated see [method_msstatsptm()].
#' @param analyses which models to return.
#' @param method label written in the `method` column.
#' @return tidy results with feature ids = PTM and PROTEIN `ProteinName`s.
#' @details Use unique `BioReplicate` labels per biological sample: MSstats treats a
#'   label shared by several conditions as the same subject (repeated measures).
#' @export
method_msstatsptm_native <- function(input, contrasts, data_type = c("LabelFree", "TMT"), MBimpute = TRUE,
                                     normalization = "equalizeMedians", moderated = FALSE,
                                     analyses = c("protein_DPA", "ptm_DPA", "ptm_DPU"),
                                     method = if (MBimpute) "msstatsptm_native" else "msstatsptm_native_noimpute") {
  data_type <- match.arg(data_type)
  if (!requireNamespace("MSstatsPTM", quietly = TRUE)) stop("Package MSstatsPTM is not installed")
  settings <- if (data_type == "TMT")
    list(method = "msstats", global_norm = !isFALSE(normalization), global_norm.PTM = !isFALSE(normalization),
         MBimpute = MBimpute, MBimpute.PTM = MBimpute)
  else list(normalization = normalization, normalization.PTM = normalization, MBimpute = MBimpute, MBimpute.PTM = MBimpute)
  gc <- msstats_fit(input, contrasts, data_type, settings, moderated)
  msstats_tidy(gc, analyses, site_map = NULL, prot_map = NULL, keep_proteins = NULL, method = method)
}

# ---- helpers ------------------------------------------------------------------------

msstats_group <- function(ann) if (!all(is.na(ann$factor))) paste(ann$condition, ann$factor, sep = ".") else ann$condition
msstats_biorep <- function(ann) if (all(is.na(ann$subject))) ann$sample else ann$subject

msstats_long_lf <- function(M, names, ann) {
  df <- data.frame(ProteinName = rep(names, times = ncol(M)), PeptideSequence = rep(names, times = ncol(M)),
                   PrecursorCharge = 2L, FragmentIon = NA_character_, ProductCharge = NA_integer_,
                   IsotopeLabelType = "L", Condition = rep(msstats_group(ann), each = nrow(M)),
                   BioReplicate = rep(msstats_biorep(ann), each = nrow(M)), Run = rep(ann$sample, each = nrow(M)),
                   Fraction = 1L, Intensity = 2^as.vector(M), stringsAsFactors = FALSE)
  df[!is.na(df$Intensity), ]
}

msstats_long_tmt <- function(M, names, ann) {
  mix <- if (is.null(ann$mixture) || all(is.na(ann$mixture))) rep("1", nrow(ann)) else ann$mixture
  # MSstatsTMT splits PSM at "_" into PeptideSequence and Charge: the peptide label must not contain "_"
  pep <- gsub("_", ".", names, fixed = TRUE)
  df <- data.frame(ProteinName = rep(names, times = ncol(M)), PeptideSequence = rep(pep, times = ncol(M)),
                   Charge = 2L, PSM = rep(paste0(pep, "_2"), times = ncol(M)),
                   Mixture = rep(mix, each = nrow(M)), TechRepMixture = 1L,
                   Run = rep(paste0(mix, "_1"), each = nrow(M)), Channel = rep(ann$sample, each = nrow(M)),
                   Condition = rep(msstats_group(ann), each = nrow(M)),
                   BioReplicate = rep(msstats_biorep(ann), each = nrow(M)),
                   Intensity = 2^as.vector(M), stringsAsFactors = FALSE)
  df[!is.na(df$Intensity), ]
}

msstats_contrast_matrix <- function(levels, contrasts) {
  cm <- t(vapply(seq_len(nrow(contrasts)), function(k) {
    v <- stats::setNames(rep(0, length(levels)), levels)
    g1 <- trimws(strsplit(contrasts$group1[k], "+", fixed = TRUE)[[1]])
    g2 <- trimws(strsplit(contrasts$group2[k], "+", fixed = TRUE)[[1]])
    if (!all(c(g1, g2) %in% levels)) stop("Contrast groups not found in MSstats GROUP levels: ",
                                         paste(setdiff(c(g1, g2), levels), collapse = ", "))
    v[g1] <- 1 / length(g1); v[g2] <- -1 / length(g2); v
  }, numeric(length(levels))))
  if (is.null(dim(cm))) cm <- matrix(cm, nrow = 1, dimnames = list(NULL, levels))
  rownames(cm) <- contrasts$name
  cm
}

# Logging goes to a temporary file. The MSstatsPTM tutorial's `use_log_file = FALSE,
# append = FALSE` is also valid; it failed in one of our first runs with
# "Assertion on '!(append & !use_log_file)' failed" for reasons we did not identify,
# and the temporary file avoids the question.
msstats_fit <- function(input, contrasts, data_type, settings, moderated) {
  logf <- tempfile("MSstatsPTM_", fileext = ".log"); on.exit(unlink(logf))
  log_args <- list(use_log_file = TRUE, append = FALSE, verbose = FALSE, log_file_path = logf)
  summ_fun <- if (data_type == "TMT") MSstatsPTM::dataSummarizationPTM_TMT else MSstatsPTM::dataSummarizationPTM
  summ <- do.call(summ_fun, c(list(input), settings, log_args))
  pl <- summ$PTM$ProteinLevelData                       # GROUP (label-free) or Condition (TMT)
  lv <- sort(unique(as.character(if ("GROUP" %in% names(pl)) pl$GROUP else pl$Condition)))
  cm <- msstats_contrast_matrix(lv, contrasts)
  log_args$append <- TRUE
  # fitted models are not kept when the installed MSstatsPTM allows it (argument added after 2.4): DiffPTM
  # only uses the comparison tables, and one model per site and protein exhausts memory on large datasets
  args <- c(list(summ, data.type = data_type, contrast.matrix = cm, moderated = moderated), log_args)
  if ("save_fitted_models" %in% names(formals(MSstatsPTM::groupComparisonPTM))) args$save_fitted_models <- FALSE
  do.call(MSstatsPTM::groupComparisonPTM, args)
}

msstats_tidy <- function(gc, analyses, site_map, prot_map, keep_proteins, method) {
  tidy <- function(df, analysis, map = NULL) {
    if (is.null(df) || !nrow(df)) return(empty_result())
    id <- as.character(df$Protein)
    if (!is.null(map)) id <- map$feature_id[match(id, map$key)]
    out <- data.frame(feature_id = id, analysis = analysis, contrast = as.character(df$Label),
                      log2FC = df$log2FC, se = df$SE, df = if (is.null(df$DF)) NA_real_ else df$DF,
                      statistic = if (is.null(df$Tvalue)) df$log2FC / df$SE else df$Tvalue, pvalue = df$pvalue,
                      adj_pvalue = df$adj.pvalue, method = method, stringsAsFactors = FALSE)
    out[!is.finite(out$log2FC), c("log2FC", "se", "df", "statistic", "pvalue", "adj_pvalue")] <- NA_real_
    if (analysis == "ptm_DPU" && "Adjusted" %in% names(df)) {
      # ADJUSTED.Model keeps sites whose protein was not found with their *unadjusted*
      # statistics (Adjusted == FALSE): drop them from DPU and re-apply BH to the rest
      na <- !df$Adjusted %in% TRUE
      out[na, c("log2FC", "se", "df", "statistic", "pvalue", "adj_pvalue")] <- NA_real_
      for (k in split(seq_len(nrow(out)), out$contrast)) {
        ok <- k[!is.na(out$pvalue[k])]
        out$adj_pvalue[ok] <- stats::p.adjust(out$pvalue[ok], "BH")
      }
    }
    out
  }
  out <- list()
  if ("ptm_DPU" %in% analyses) out$dpu <- tidy(gc$ADJUSTED.Model, "ptm_DPU", site_map)
  if ("ptm_DPA" %in% analyses) out$dpa <- tidy(gc$PTM.Model, "ptm_DPA", site_map)
  if ("protein_DPA" %in% analyses) {
    pr <- tidy(gc$PROTEIN.Model, "protein_DPA", prot_map)
    if (!is.null(keep_proteins)) pr <- pr[pr$feature_id %in% keep_proteins, ]
    out$prot <- pr
  }
  res <- do.call(rbind, out); rownames(res) <- NULL
  res[!is.na(res$feature_id), ]
}
