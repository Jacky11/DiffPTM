#' Match PTM features to their parent protein in the total proteome
#'
#' Two strategies, usable alone or combined (`"hybrid"`: sequence first, then ID for
#' features without a sequence hit):
#'
#' * `"id"`: a PTM feature matches every protein-group row of the proteome that shares at
#'   least one accession with `meta$protein_ids` (members split on the assay separator;
#'   isoform suffixes `-2` are optionally stripped).
#' * `"sequence"`: the stripped peptide sequence is searched in the protein sequences of a
#'   FASTA file (Aho-Corasick when `AhoCorasickTrie` is installed, fixed `grepl` loop
#'   otherwise); the accession(s) found are then mapped to proteome rows by ID. This is the
#'   strategy recommended by Wolski et al. (submitted to Methods Mol Biol, §1.6) because it bypasses differences
#'   in protein inference between the PTM and total-proteome searches.
#'
#' When a PTM feature maps to several quantified protein rows, the parent-protein
#' intensity is the per-sample mean of their log2 intensities (`shared = "mean"`), the
#' first row (`"first"`) or the feature is dropped (`"drop"`).
#'
#' The result is stored in `x$ptm$meta$matched_protein` (feature_id(s) of the proteome
#' rows, separated by `;`) and `x$ptm$protein_intensity`, a matrix aligned on the PTM
#' features (rows) and the **protein** samples (columns). Run matching *after*
#' normalisation / imputation of the protein assay so that the copied values are final.
#'
#' @param x a `ptm_experiment` with both assays.
#' @param method `"id"`, `"sequence"` or `"hybrid"`.
#' @param fasta path to a FASTA file (required for `"sequence"` / `"hybrid"`).
#' @param strip_isoforms remove `-n` suffixes from accessions before matching.
#' @param shared handling of PTM features matched to several protein rows.
#' @param drop_unmatched remove PTM features without parent protein (default keeps them:
#'   they still get DPA results, and `NA` for DPU).
#' @export
match_ptm_to_protein <- function(x, method = c("id", "sequence", "hybrid"), fasta = NULL,
                                 strip_isoforms = TRUE, shared = c("mean", "first", "drop"),
                                 drop_unmatched = FALSE) {
  method <- match.arg(method); shared <- match.arg(shared)
  if (is.null(x$ptm) || is.null(x$protein)) stop("Both assays are required for matching")

  prot_meta <- x$protein$meta
  prot_members <- split_ids(prot_meta$protein_ids, x$protein$sep, strip_isoforms)
  lookup <- data.frame(acc = unlist(prot_members),
                       row = rep(seq_along(prot_members), lengths(prot_members)),
                       stringsAsFactors = FALSE)
  lookup <- lookup[!duplicated(lookup$acc), ]
  acc_to_row <- stats::setNames(lookup$row, lookup$acc)

  ptm_members <- split_ids(x$ptm$meta$protein_ids, x$ptm$sep, strip_isoforms)
  rows_by_id <- lapply(ptm_members, function(acc) unique(unname(acc_to_row[acc[acc %in% names(acc_to_row)]])))

  rows_by_seq <- NULL
  if (method %in% c("sequence", "hybrid")) {
    if (is.null(fasta)) stop("`fasta` is required for sequence-based matching")
    seqs <- read_fasta(fasta, strip_isoforms = strip_isoforms)
    pep <- toupper(x$ptm$meta$stripped_sequence)
    hits <- find_peptides(pep, seqs)       # list: per peptide, accessions containing it
    rows_by_seq <- lapply(hits, function(acc) unique(unname(acc_to_row[acc[acc %in% names(acc_to_row)]])))
  }
  rows <- switch(method,
    id = rows_by_id,
    sequence = rows_by_seq,
    hybrid = mapply(function(s, i) if (length(s)) s else i, rows_by_seq, rows_by_id, SIMPLIFY = FALSE)
  )
  n_match <- lengths(rows)
  if (shared == "first") rows <- lapply(rows, function(r) r[seq_len(min(1, length(r)))])
  keep <- if (shared == "drop") n_match <= 1 else rep(TRUE, length(n_match))
  if (drop_unmatched) keep <- keep & n_match >= 1

  P <- x$protein$intensity
  mp <- t(vapply(rows, function(r) {
    if (!length(r)) return(rep(NA_real_, ncol(P)))
    if (length(r) == 1) return(P[r, ])
    colMeans(P[r, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(P))))
  mp[is.nan(mp)] <- NA_real_
  colnames(mp) <- colnames(P); rownames(mp) <- x$ptm$meta$feature_id

  x$ptm$meta$matched_protein <- vapply(rows, function(r) paste(prot_meta$feature_id[r], collapse = ";"), "")
  x$ptm$meta$n_matched_protein <- n_match
  # stored inside the assay so that later subsetting keeps it aligned
  x$ptm$protein_intensity <- mp
  if (!is.null(x$protein$n_features)) {                # precision weights travel with the parent values
    NF <- x$protein$n_features
    x$ptm$protein_n_features <- t(vapply(rows, function(r) if (!length(r)) rep(NA_real_, ncol(NF)) else colSums(NF[r, , drop = FALSE]), numeric(ncol(NF))))
    dimnames(x$ptm$protein_n_features) <- dimnames(mp)
  }
  x$ptm <- subset_assay(x$ptm, keep)
  add_log(x, sprintf("Matched PTM features to parent proteins by %s (%s); %d/%d features matched, shared peptides: %s%s.",
                     method, if (strip_isoforms) "isoform suffixes ignored" else "exact accessions",
                     sum(n_match >= 1), length(n_match), shared,
                     if (drop_unmatched) "; unmatched features removed" else ""))
}

split_ids <- function(ids, sep, strip_isoforms = TRUE) {
  ids[is.na(ids)] <- ""
  out <- strsplit(ids, sep, fixed = TRUE)
  out <- lapply(out, function(v) { v <- trimws(v); v <- sub("^(?:sp|tr)\\|([^|]+)\\|.*$", "\\1", v, perl = TRUE); v[nzchar(v)] })
  if (strip_isoforms) out <- lapply(out, function(v) sub("-\\d+$", "", v))
  out
}

#' Read a FASTA file into a named character vector (accession -> sequence)
#' @keywords internal
read_fasta <- function(path, strip_isoforms = TRUE) {
  if (requireNamespace("Biostrings", quietly = TRUE)) {
    aa <- Biostrings::readAAStringSet(path)
    seqs <- as.character(aa); nm <- names(aa)
  } else {
    lines <- readLines(path)
    hdr <- grepl("^>", lines)
    idx <- cumsum(hdr)
    nm <- sub("^>", "", lines[hdr])
    seqs <- vapply(split(lines[!hdr], idx[!hdr]), paste, "", collapse = "")
  }
  acc <- sub("^(?:sp|tr)\\|([^|]+)\\|.*$", "\\1", nm, perl = TRUE)
  acc <- sub(" .*$", "", acc)
  if (strip_isoforms) acc <- sub("-\\d+$", "", acc)
  seqs <- toupper(gsub("[^A-Z]", "", toupper(seqs)))
  stats::setNames(seqs, acc)
}

#' Find which protein sequences contain each peptide
#' @return list (one element per peptide) of accessions.
#' @keywords internal
find_peptides <- function(peptides, seqs) {
  res <- vector("list", length(peptides))
  ok <- !is.na(peptides) & nzchar(peptides)
  uniq <- unique(peptides[ok])
  if (!length(uniq)) return(res)
  if (requireNamespace("AhoCorasickTrie", quietly = TRUE)) {
    # I/L are indistinguishable by MS: search on an I->L normalised alphabet
    norm <- function(s) chartr("I", "L", s)
    texts <- stats::setNames(norm(unname(seqs)), names(seqs))
    hits <- AhoCorasickTrie::AhoCorasickSearch(norm(uniq), texts, alphabet = "aminoacid", groupByKeyword = FALSE)
    # hits: one element per text; each element is a list of matches (list(Keyword, Offset)) or a
    # data.frame, depending on the package version. Texts are identified by name when the
    # package keeps them, by position otherwise.
    per_pep <- parse_ac_hits(hits, norm(uniq), names(seqs))
    found <- if (is.null(per_pep)) lapply(norm(peptides[ok]), function(p) names(seqs)[grepl(p, texts, fixed = TRUE)])
             else per_pep[norm(peptides[ok])]
  } else {
    nseq <- chartr("I", "L", seqs)
    found <- lapply(chartr("I", "L", peptides[ok]), function(p) names(seqs)[grepl(p, nseq, fixed = TRUE)])
  }
  res[ok] <- lapply(found, function(v) unique(v %||% character(0)))
  res
}

#' Turn AhoCorasickTrie::AhoCorasickSearch output into a keyword -> text-names list
#'
#' `hits` has one element per text, each a list of matches (`list(Keyword, Offset)`) or a
#' data.frame with a `Keyword` column. Returns `NULL` when texts cannot be identified
#' (caller then falls back to `grepl`).
#' @keywords internal
parse_ac_hits <- function(hits, keywords, text_names) {
  nm <- if (!is.null(names(hits)) && all(nzchar(names(hits)))) names(hits) else if (length(hits) == length(text_names)) text_names else NULL
  if (is.null(nm)) return(NULL)
  per_kw <- stats::setNames(vector("list", length(keywords)), keywords)
  for (i in seq_along(hits)) {
    h <- hits[[i]]
    kws <- if (is.data.frame(h)) h$Keyword else if (is.list(h)) unlist(lapply(h, function(z) z[["Keyword"]])) else NULL
    for (k in unique(kws)) per_kw[[k]] <- c(per_kw[[k]], nm[i])
  }
  per_kw
}
