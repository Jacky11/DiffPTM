#' DiffPTM experiment object
#'
#' A `ptm_experiment` holds two assays (`protein`: total proteome, `ptm`: modified
#' peptidoforms or sites), one sample annotation per assay, processing parameters and
#' a processing log used to generate the Materials & Methods paragraph.
#'
#' @param protein,ptm assays created with [new_assay()] (either may be `NULL` while importing).
#' @param annotation named list with elements `protein` and/or `ptm`, each a data.frame
#'   created with [make_annotation()].
#' @param params list of user parameters (kept for reporting).
#' @return an object of class `ptm_experiment`.
#' @export
new_ptm_experiment <- function(protein = NULL, ptm = NULL, annotation = list(), params = list()) {
  x <- structure(
    list(protein = protein, ptm = ptm, annotation = annotation,
         params = params, log = character(0), results = list()),
    class = "ptm_experiment"
  )
  validate_ptm_experiment(x)
}

#' Create an assay
#'
#' @param intensity numeric matrix, features in rows, samples in columns. Values are
#'   expected on the log2 scale; set `log2_transform = TRUE` to transform raw intensities
#'   (zeros are converted to `NA` first, as in the original app).
#' @param meta data.frame with one row per feature. Must contain `feature_id` (unique).
#'   For PTM assays the following columns are used when present: `protein_ids`
#'   (members separated by `sep`), `stripped_sequence`, `modified_sequence`, `site_key`,
#'   `residue`, `position`, `multiplicity`, `localization_prob`, `ptm_type`, `gene`.
#' @param log2_transform logical.
#' @param sep separator between protein identifiers in `meta$protein_ids`.
#' @return a list of class `ptm_assay` with elements `intensity`, `raw_detected`
#'   (logical matrix of observed values before any imputation), `meta`, `sep`.
#' @export
new_assay <- function(intensity, meta, log2_transform = FALSE, sep = ";") {
  intensity <- as.matrix(intensity)
  storage.mode(intensity) <- "double"
  if (log2_transform) {
    intensity[!is.na(intensity) & intensity <= 0] <- NA_real_
    intensity <- log2(intensity)
  }
  if (!"feature_id" %in% names(meta)) stop("`meta` must contain a `feature_id` column")
  if (anyDuplicated(meta$feature_id)) stop("`meta$feature_id` must be unique")
  if (nrow(meta) != nrow(intensity)) stop("`meta` and `intensity` must have the same number of rows")
  rownames(intensity) <- meta$feature_id
  structure(
    list(intensity = intensity, raw_detected = !is.na(intensity),
         meta = as.data.frame(meta, stringsAsFactors = FALSE), sep = sep),
    class = "ptm_assay"
  )
}

validate_ptm_experiment <- function(x) {
  for (a in c("protein", "ptm")) {
    if (is.null(x[[a]])) next
    if (!inherits(x[[a]], "ptm_assay")) stop(sprintf("`%s` must be created with new_assay()", a))
    ann <- x$annotation[[a]]
    if (!is.null(ann)) {
      missing <- setdiff(ann$sample, colnames(x[[a]]$intensity))
      if (length(missing)) stop(sprintf("Annotation for `%s` refers to unknown samples: %s",
                                        a, paste(missing, collapse = ", ")))
      # keep only annotated samples, in annotation order
      x[[a]]$intensity    <- x[[a]]$intensity[, ann$sample, drop = FALSE]
      x[[a]]$raw_detected <- x[[a]]$raw_detected[, ann$sample, drop = FALSE]
    }
  }
  x
}

#' @export
print.ptm_experiment <- function(x, ...) {
  cli::cli_h1("ptm_experiment")
  for (a in c("protein", "ptm")) {
    if (is.null(x[[a]])) { cli::cli_text("{a}: <none>"); next }
    d <- dim(x[[a]]$intensity)
    cli::cli_text("{a}: {d[1]} features x {d[2]} samples; ",
                  "{round(100 * mean(is.na(x[[a]]$intensity)), 1)}% missing")
    ann <- x$annotation[[a]]
    if (!is.null(ann)) cli::cli_text("   groups: {paste(unique(groups_of(x, a)), collapse = ', ')}")
  }
  if (length(x$log)) { cli::cli_h2("log"); cli::cli_ol(x$log) }
  invisible(x)
}

#' Append a processing step to the log
#' @keywords internal
add_log <- function(x, msg) {
  x$log <- c(x$log, msg)
  x
}

#' Group label (condition, or condition:factor) of each sample of an assay
#' @keywords internal
groups_of <- function(x, assay) {
  ann <- x$annotation[[assay]]
  if (is.null(ann)) stop(sprintf("No annotation for assay `%s`", assay))
  if ("factor" %in% names(ann) && !all(is.na(ann$factor))) {
    factor(paste(ann$condition, ann$factor, sep = "."),
           levels = unique(paste(ann$condition, ann$factor, sep = ".")))
  } else {
    factor(ann$condition, levels = unique(ann$condition))
  }
}

#' Sample annotation
#'
#' @param sample column names of the intensity matrix.
#' @param condition biological condition of each sample.
#' @param factor optional second factor (time point, batch...). `NULL` or `NA` = none.
#' @param subject optional subject / pairing identifier used as a blocking factor.
#' @param mixture optional TMT plex (mixture) of each sample, used as a fixed batch effect.
#' @return data.frame with columns `sample, condition, factor, subject, mixture`.
#' @export
make_annotation <- function(sample, condition, factor = NULL, subject = NULL, mixture = NULL) {
  n <- length(sample)
  if (length(condition) != n) stop("`condition` must have one value per sample")
  if (is.null(factor))  factor  <- rep(NA_character_, n)
  if (is.null(subject)) subject <- rep(NA_character_, n)
  if (is.null(mixture)) mixture <- rep(NA_character_, n)
  ann <- data.frame(sample = as.character(sample), condition = as.character(condition),
                    factor = as.character(factor), subject = as.character(subject),
                    mixture = as.character(mixture), stringsAsFactors = FALSE)
  if (anyDuplicated(ann$sample)) stop("duplicated sample names in annotation")
  tab <- table(ann$condition, ann$factor, useNA = "ifany")
  if (any(tab[tab > 0] < 2)) warning("Some groups have fewer than 2 replicates")
  ann
}

#' Guess an annotation from sample names
#'
#' Sample names such as `WT_t0_rep1` are split with `pattern` into tokens; the user then
#' maps token positions to condition / factor / subject.
#'
#' @param sample character vector.
#' @param pattern regular expression separating tokens.
#' @param condition_pos,factor_pos,subject_pos token index (1-based), a vector of indices, a
#'   string such as `"1-2"` (tokens joined with `_`), `-1` (all tokens but the last, i.e. the
#'   sample name without its replicate number), `0` or `NULL`.
#' @export
guess_annotation <- function(sample, pattern = "[_ .-]", condition_pos = 1,
                             factor_pos = NULL, subject_pos = NULL) {
  tok <- strsplit(sample, pattern)
  pick <- function(pos) {
    pos <- parse_token_positions(pos)
    if (is.null(pos)) return(NULL)
    vapply(tok, function(t) {
      idx <- if (all(pos < 0)) seq_len(max(length(t) + min(pos), 0)) else pos[pos <= length(t)]   # -1: all but the last token
      v <- t[idx]; if (length(v)) paste(v, collapse = "_") else NA_character_ }, "")
  }
  make_annotation(sample, pick(condition_pos), pick(factor_pos), pick(subject_pos))
}

#' Token positions given as a number, a vector or a string such as `"1-2"` or `"1,3"`
#' @keywords internal
parse_token_positions <- function(pos) {
  if (is.null(pos) || !length(pos)) return(NULL)
  if (is.character(pos)) {
    pos <- trimws(pos); if (!nzchar(pos) || pos == "0") return(NULL)
    pos <- unlist(lapply(strsplit(pos, ",", fixed = TRUE)[[1]], function(r) {
      r <- trimws(r)
      if (grepl("^-\\d+$", r)) return(as.integer(r))                        # negative: all tokens but the last n
      r <- as.integer(strsplit(r, "-", fixed = TRUE)[[1]]); if (length(r) == 2) seq(r[1], r[2]) else r }))
  }
  pos <- as.integer(pos); pos <- pos[!is.na(pos) & pos != 0]
  if (!length(pos)) NULL else pos
}

#' Rename or select the samples of an assay
#'
#' `rename_samples()` selects and renames samples with a regular expression whose first
#' capture group is the new name: `"^Phospho_(.+)$"` keeps the samples `Phospho_mock_24h_1`,
#' ... and names them `mock_24h_1`, ..., so that the PTM and proteome samples of a combined
#' MaxQuant search (whose site table also lists the proteome runs) get the same names and
#' are paired. Samples that do not match the pattern are dropped (`drop_unmatched = TRUE`)
#' or kept under their old name. `keep_samples()` restricts an assay to the given samples.
#'
#' @param a a `ptm_assay`.
#' @param pattern regular expression with one capture group.
#' @param drop_unmatched drop the samples that do not match `pattern`.
#' @param samples sample names to keep.
#' @export
rename_samples <- function(a, pattern, drop_unmatched = TRUE) {
  if (is.null(pattern) || !nzchar(pattern)) return(a)
  old <- colnames(a$intensity)
  hit <- grepl(pattern, old, perl = TRUE)
  if (!any(hit)) stop("No sample name matches the pattern")
  if (drop_unmatched) a <- keep_samples(a, old[hit])
  old <- colnames(a$intensity); hit <- grepl(pattern, old, perl = TRUE)
  new <- old; new[hit] <- sub(pattern, "\\1", old[hit], perl = TRUE)
  if (anyDuplicated(new)) stop("Renaming gives duplicated sample names: ", paste(unique(new[duplicated(new)]), collapse = ", "))
  colnames(a$intensity) <- colnames(a$raw_detected) <- new
  if (!is.null(a$n_features)) colnames(a$n_features) <- new
  a
}

#' @rdname rename_samples
#' @export
keep_samples <- function(a, samples) {
  samples <- intersect(samples, colnames(a$intensity))
  if (!length(samples)) stop("No sample left")
  a$intensity <- a$intensity[, samples, drop = FALSE]; a$raw_detected <- a$raw_detected[, samples, drop = FALSE]
  if (!is.null(a$n_features)) a$n_features <- a$n_features[, samples, drop = FALSE]
  a
}

#' Read an annotation table (tsv/csv)
#'
#' Required columns: `sample`, `condition`; optional: `factor`, `subject`, `mixture`.
#' @param path file path.
#' @return see [make_annotation()].
#' @export
read_annotation <- function(path) {
  d <- data.table::fread(path, data.table = FALSE, colClasses = "character")
  names(d) <- tolower(names(d))
  need <- c("sample", "condition")
  if (!all(need %in% names(d))) stop("annotation file needs columns: sample, condition")
  make_annotation(d$sample, d$condition, d[["factor"]], d[["subject"]], d[["mixture"]])
}
