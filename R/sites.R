#' Collapse peptidoforms to modification sites
#'
#' Features sharing the same `site_key` (default: first protein accession + residue +
#' position; for MaxQuant sites tables `positions_within_proteins` of the leading protein)
#' are aggregated per sample. Aggregation is done on the linear scale for `"sum"` and on
#' the log2 scale for `"max"` / `"median"`. Metadata of the first member is kept and the
#' number of collapsed peptidoforms is stored in `n_peptidoforms`.
#'
#' @param x a `ptm_experiment`.
#' @param method `"sum"`, `"max"` or `"median"`.
#' @param site_key optional character vector overriding the computed key.
#' @export
collapse_to_sites <- function(x, method = c("sum", "max", "median"), site_key = NULL) {
  method <- match.arg(method)
  if (is.null(x$ptm)) return(x)
  meta <- x$ptm$meta
  key <- site_key %||% make_site_key(meta, x$ptm$sep)
  if (all(is.na(key))) stop("Cannot build site keys: residue/position information missing")
  key[is.na(key)] <- meta$feature_id[is.na(key)]        # unresolvable -> keep as is
  M <- x$ptm$intensity
  idx <- split(seq_len(nrow(M)), factor(key, levels = unique(key)))
  agg <- t(vapply(idx, function(i) {
    sub <- M[i, , drop = FALSE]
    switch(method,
      sum    = { v <- log2(colSums(2^sub, na.rm = TRUE)); v[colSums(!is.na(sub)) == 0] <- NA; v },
      max    = apply(sub, 2, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)),
      median = apply(sub, 2, function(v) if (all(is.na(v))) NA_real_ else stats::median(v, na.rm = TRUE)))
  }, numeric(ncol(M))))
  colnames(agg) <- colnames(M)
  new_meta <- meta[vapply(idx, `[`, 1L, 1L), , drop = FALSE]
  new_meta$n_peptidoforms <- lengths(idx)
  new_meta$site_key <- names(idx)
  new_meta$feature_id <- make.unique(names(idx))
  x$ptm <- new_assay(agg, new_meta, log2_transform = FALSE, sep = x$ptm$sep)
  add_log(x, sprintf("Collapsed %d peptidoforms to %d sites (%s).", nrow(M), nrow(agg), method))
}

make_site_key <- function(meta, sep = ";") {
  lead <- vapply(split_ids(meta$protein_ids, sep, strip_isoforms = FALSE), function(v) if (length(v)) v[1] else NA_character_, "")
  pos <- meta$position
  if (all(is.na(pos)) && "positions_within_proteins" %in% names(meta))
    pos <- vapply(strsplit(as.character(meta$positions_within_proteins), sep, fixed = TRUE), `[`, "", 1)
  res <- meta$residue
  k <- paste0(lead, "_", ifelse(is.na(res), "", res), pos)
  k[is.na(lead) | is.na(pos)] <- NA_character_
  if ("ptm_type" %in% names(meta) && !all(is.na(meta$ptm_type))) k <- paste0(k, "[", meta$ptm_type, "]")
  k
}

#' Group peptidoforms into FragPipe-style "multisite features"
#'
#' All peptidoforms carrying the same *set* of modification sites on the same protein get
#' the same `multisite_key`, regardless of missed cleavages or additional variable
#' modifications. The key is added to `meta`; call [collapse_to_sites()] with
#' `site_key = x$ptm$meta$multisite_key` to aggregate.
#'
#' @param x a `ptm_experiment`.
#' @param position_sep separator between positions when a feature carries several
#'   (MaxQuant "Positions within proteins" uses `;` for proteins and the sites table has one
#'   row per site; Spectronaut `PTM.CollapseKey` encodes `P12345_S123_M1`).
#' @export
multisite_features <- function(x, position_sep = ";") {
  if (is.null(x$ptm)) return(x)
  meta <- x$ptm$meta
  lead <- vapply(split_ids(meta$protein_ids, x$ptm$sep, FALSE), function(v) if (length(v)) v[1] else NA_character_, "")
  sites <- lapply(strsplit(as.character(meta$position), position_sep, fixed = TRUE),
                  function(v) paste(sort(unique(trimws(v))), collapse = "+"))
  meta$multisite_key <- paste0(lead, "_", unlist(sites))
  x$ptm$meta <- meta
  add_log(x, sprintf("Built %d multisite features from %d peptidoforms.", length(unique(meta$multisite_key)), nrow(meta)))
}

#' Merge PTM rows that report the same evidence for several proteins
#'
#' Site reports of some software (e.g. Spectronaut PTM site reports) contain one row per
#' protein accession when a modified peptide maps to several proteins or isoforms. These
#' rows carry identical quantities, so keeping them all counts the same measurement
#' several times: site numbers are inflated and the multiple-testing correction is
#' distorted. Rows sharing `key` (default: `modified_sequence`, plus `multiplicity` when
#' present) **and** identical intensities in every sample are merged into one feature.
#' The kept row is the one whose accession looks reviewed (UniProt 6-character
#' accession), otherwise the first; `protein_ids` becomes the union of the accessions so
#' that [match_ptm_to_protein()] still finds the parent protein(s). Rows with the same
#' key but different quantities are left untouched.
#'
#' @param x a `ptm_experiment`.
#' @param key optional character vector (one value per PTM row) overriding the default key.
#' @param tol tolerance on log2 intensities for rows to be considered identical.
#' @export
merge_redundant_ptm <- function(x, key = NULL, tol = 1e-8) {
  if (is.null(x$ptm)) return(x)
  meta <- x$ptm$meta; M <- x$ptm$intensity
  if (is.null(key)) {
    if (!"modified_sequence" %in% names(meta) || all(is.na(meta$modified_sequence)))
      stop("`modified_sequence` is missing: supply `key`")
    key <- as.character(meta$modified_sequence)
    if ("multiplicity" %in% names(meta)) key <- paste(key, meta$multiplicity, sep = "|")
  }
  qkey <- apply(round(M / tol) * tol, 1, function(v) paste(ifelse(is.na(v), "NA", format(v, digits = 15)), collapse = ","))
  grp <- paste(key, qkey, sep = "#")
  grp[is.na(key)] <- paste0(".row", which(is.na(key)))
  idx <- split(seq_len(nrow(M)), factor(grp, levels = unique(grp)))
  if (all(lengths(idx) == 1)) return(add_log(x, "No redundant PTM rows found."))
  first_acc <- vapply(split_ids(meta$protein_ids, x$ptm$sep, strip_isoforms = FALSE),
                      function(v) if (length(v)) v[1] else NA_character_, "")
  reviewed <- grepl("^([OPQ][0-9][A-Z0-9]{3}[0-9])(-\\d+)?$", first_acc)
  keep <- vapply(idx, function(i) if (any(reviewed[i])) i[reviewed[i]][1] else i[1], 1L)
  acc <- vapply(idx, function(i) paste(unique(unlist(split_ids(meta$protein_ids[i], x$ptm$sep, FALSE))), collapse = x$ptm$sep), "")
  new_meta <- meta[keep, , drop = FALSE]
  new_meta$protein_ids <- acc
  new_meta$n_redundant <- lengths(idx)
  x$ptm <- new_assay(M[keep, , drop = FALSE], new_meta, log2_transform = FALSE, sep = x$ptm$sep)
  add_log(x, sprintf("Merged %d PTM rows reporting the same peptide-level evidence for several proteins into %d features.",
                     nrow(M), length(idx)))
}

#' Derive site residues and protein positions for peptidoform tables
#'
#' Precursor or modified-peptide tables (DIA-NN `report.pr_matrix.tsv`, FragPipe
#' `combined_modified_peptide.tsv`, generic tables) give the modified sequence but not the
#' position of the site in the protein. This function reads the position of each
#' modification inside the peptide from the modified sequence (a residue followed by an
#' annotation such as `(UniMod:21)`, `[79.9663]`, `(ph)` or `(Phospho (STY))`) and, when a
#' FASTA file is given, locates the peptide in its protein (first accession of
#' `protein_ids`, I/L treated as equivalent) to obtain the position in the protein.
#' Residues annotated with fixed or label modifications (UniMod:4 carbamidomethyl, UniMod:35
#' oxidation, UniMod:1 acetyl, 57.02, 15.99, 42.01, TMT 229.16 / 304.21) are ignored.
#' Multiply modified peptides get `;`-separated residues and positions, and a
#' `multiplicity` equal to the number of sites.
#'
#' @param x a `ptm_experiment`.
#' @param fasta path to the FASTA of the search (optional; without it only the residue and
#'   the position within the peptide are filled).
#' @param ignore regex of annotations to ignore.
#' @export
locate_sites <- function(x, fasta = NULL,
                         ignore = "UniMod:4\\b|UniMod:35\\b|UniMod:1\\b|UniMod:737|UniMod:2016|57\\.02|15\\.99|42\\.01|229\\.16|304\\.20|304\\.21|Carbamidomethyl|Oxidation|Acetyl|TMT") {
  if (is.null(x$ptm)) return(x)
  meta <- x$ptm$meta
  ms <- as.character(meta$modified_sequence); ss <- as.character(meta$stripped_sequence)
  todo <- !is.na(ms) & (is.na(meta$position) | !nzchar(meta$position))
  if (!any(todo)) return(add_log(x, "Site positions already present; nothing located."))
  parsed <- lapply(ms[todo], function(m) {
    m <- sub("^[A-Z_-]?\\.", "", m); m <- sub("\\.[A-Z_-]?$", "", m)        # PD-style flanking "K.PEPTIDE.R"
    m <- gsub("^\\(.*?\\)|^\\[.*?\\]|^_", "", m)                            # N-terminal annotation
    toks <- regmatches(m, gregexpr("[A-Za-z](?:\\([^()]*(?:\\([^()]*\\))?[^()]*\\)|\\[[^]]*\\])?", m, perl = TRUE))[[1]]
    toks <- toks[grepl("^[A-Z]", toks)]
    ann <- sub("^[A-Z]", "", toks); aa <- substr(toks, 1, 1)
    is_site <- nzchar(ann) & !grepl(ignore, ann, perl = TRUE)
    list(stripped = paste(aa, collapse = ""), idx = which(is_site), residue = aa[is_site])
  })
  stripped <- vapply(parsed, `[[`, "", "stripped")
  ss[todo & (is.na(ss) | !nzchar(ss))] <- stripped[(is.na(ss) | !nzchar(ss))[todo]]
  meta$stripped_sequence <- ss
  idx <- lapply(parsed, `[[`, "idx"); resid <- lapply(parsed, `[[`, "residue")
  pep_pos <- vapply(idx, function(i) paste(i, collapse = ";"), "")
  offset <- rep(NA_integer_, sum(todo)); located <- 0L
  if (!is.null(fasta)) {
    seqs <- read_fasta(fasta)
    first <- vapply(split_ids(meta$protein_ids[todo], x$ptm$sep), function(v) if (length(v)) v[1] else NA_character_, "")
    for (k in seq_along(first)) {
      s <- seqs[first[k]]; p <- ss[todo][k]
      if (is.na(s) || is.na(p) || !nzchar(p)) next
      hit <- regexpr(chartr("I", "L", p), chartr("I", "L", s), fixed = TRUE)
      if (hit > 0) { offset[k] <- hit - 1L; located <- located + 1L }
    }
  }
  pos_out <- mapply(function(i, off) if (!length(i)) NA_character_ else if (is.na(off)) NA_character_ else paste(i + off, collapse = ";"), idx, offset)
  res_out <- vapply(resid, function(r) if (length(r)) paste(r, collapse = ";") else NA_character_, "")
  meta$residue[todo] <- res_out
  meta$position[todo] <- if (is.null(fasta)) NA_character_ else pos_out
  meta$position_in_peptide <- NA_character_; meta$position_in_peptide[todo] <- ifelse(nzchar(pep_pos), pep_pos, NA)
  n_sites <- vapply(idx, length, 1L)
  if (all(is.na(meta$multiplicity))) meta$multiplicity[todo] <- ifelse(n_sites > 0, as.character(n_sites), NA)
  meta$label <- site_label(meta, sep = x$ptm$sep)
  x$ptm$meta <- meta
  add_log(x, sprintf("Located modification sites in %d peptidoforms from the modified sequence%s.", sum(todo),
                     if (is.null(fasta)) " (residue only; give the FASTA for protein positions)" else sprintf("; %d peptides found in the FASTA", located)))
}
