#' Suggest the experimental design from the imported tables
#'
#' Proposes, without any input from the user, what the Design step of the interface needs. The
#' rules use the data and the structure of the run names, not particular words such as
#' "Phospho", so that they apply to MaxQuant, Spectronaut, DIA-NN and other exports:
#'
#' 1. **Runs without data** are excluded (e.g. the columns of a combined MaxQuant search for
#'    runs of another enrichment, which are empty in a site table).
#' 2. **Combined searches.** When both tables list the same runs (one search of the enriched and
#'    the total-proteome runs, as MaxQuant or DIA-NN can do), the PTM-enriched runs are those
#'    that quantify many more sites than the others (a gap of more than three-fold in the number
#'    of quantified sites); the other runs are excluded from the PTM table and the enriched runs
#'    from the proteome table.
#' 3. **Pairing.** A PTM run and a proteome run are the same biological sample when their names
#'    are the same once the parts that only mark the fraction are removed: parts found in every
#'    run name of one table and in few of the other (`Phospho`, `P`, `Proteome`, `T`, `PTM`,
#'    `FP`, ...). If exact matches are not found, the most similar names are paired when the
#'    match is unambiguous (e.g. a run re-injected as `..._R2_Octobre`). Paired runs receive the
#'    same sample name, which is what the usage test uses to pair them.
#' 4. **Conditions** are read from the sample names: the part that varies like a replicate
#'    number (`1`, `R1`, `rep2`, ...) is removed, as is a leading run number (`01_NS_R1`), and
#'    the rest is the condition. The position of the replicate part is chosen so that every
#'    condition has at least two samples (`R1-A` -> `A`, `mock_24h_1` -> `mock_24h`).
#'
#' @param protein,ptm `ptm_assay` objects (from [import_protein()] and [import_ptm()]).
#' @param pattern regular expression separating the parts of a run name.
#' @return a list with `annotation`, one data frame per table with columns `run` (column name in
#'   the imported table), `sample` (name used in the analysis; `""` for an excluded run),
#'   `condition` (`""` for an excluded run), `factor`, `subject` and `mixture`; and `messages`,
#'   a plain-language summary of what was done.
#' @examples
#' sim <- simulate_ptm_features(n_rep = 3, n_proteins = 100, n_sites = 80, seed = 1)
#' x <- sim_to_experiment(sim)
#' ptm <- set_run_names(x$ptm, paste0("P_", colnames(x$ptm$intensity)))
#' prot <- set_run_names(x$protein, paste0("T_", colnames(x$protein$intensity)))
#' s <- suggest_design(prot, ptm)
#' s$annotation$ptm
#' s$messages
#' @export
suggest_design <- function(protein, ptm, pattern = "[_ .-]") {
  runs <- list(protein = colnames(protein$intensity), ptm = colnames(ptm$intensity))
  nq <- list(protein = colSums(!is.na(protein$intensity)), ptm = colSums(!is.na(ptm$intensity)))
  lab <- c(protein = "proteome", ptm = "PTM")
  msg <- character()
  short <- function(v, n = 3) paste0(paste(utils::head(v, n), collapse = ", "), if (length(v) > n) ", ..." else "")
  use <- lapply(nq, function(n) n > 0)
  for (k in names(use)) {
    if (!any(use[[k]])) { use[[k]][] <- TRUE; next }
    if (any(!use[[k]])) msg <- c(msg, sprintf("%s table: %d of %d runs have no quantified value and are excluded (%s).",
                                              lab[[k]], sum(!use[[k]]), length(use[[k]]), short(runs[[k]][!use[[k]]])))
  }
  # combined search: the same runs in both tables
  common <- intersect(runs$protein[use$protein], runs$ptm[use$ptm])
  if (length(common) >= 2) {
    cnt <- nq$ptm[use$ptm]
    high <- gap_split(cnt)
    if (!is.null(high)) {
      low <- names(cnt)[!high & names(cnt) %in% runs$protein]       # few sites, and present in the proteome table
      if (length(low)) {
        use$ptm[runs$ptm %in% low] <- FALSE
        msg <- c(msg, sprintf("Combined search: %d runs quantify many more sites than the others and are taken as the PTM-enriched runs (%s); the %d others (%s) are taken as total-proteome runs.",
                              sum(high), short(names(cnt)[high]), length(low), short(low)))
        enr <- names(cnt)[high]
        if (sum(use$protein & !runs$protein %in% enr) >= 2) use$protein[runs$protein %in% enr] <- FALSE
      }
    }
  }
  # pairing and sample names
  pr <- runs$protein[use$protein]; pt <- runs$ptm[use$ptm]
  pm <- pair_runs(pr, pt, pattern)
  if (length(pm$removed$protein) || length(pm$removed$ptm)) {
    rm <- c(if (length(pm$removed$protein)) sprintf("'%s' in the proteome run names", paste(pm$removed$protein, collapse = "', '")),
            if (length(pm$removed$ptm)) sprintf("'%s' in the PTM run names", paste(pm$removed$ptm, collapse = "', '")))
    msg <- c(msg, sprintf("Sample names: %s removed, as these parts only mark the fraction.", paste(rm, collapse = " and ")))
  }
  if (pm$n_fuzzy) msg <- c(msg, sprintf("%d PTM runs paired with the most similar proteome run name (e.g. %s with %s).",
                                        pm$n_fuzzy, pm$fuzzy_example[1], pm$fuzzy_example[2]))
  paired <- intersect(pm$sample$protein, pm$sample$ptm)
  lone <- pt[!pm$sample$ptm %in% paired]
  msg <- c(msg, if (length(paired))
    sprintf("%d of %d PTM samples are paired with a proteome sample (same sample name; used by the usage test)%s.",
            length(paired), length(pt), if (length(lone)) sprintf("; not paired: %s", short(lone)) else "")
    else "No PTM run could be paired with a proteome run: usage will be tested without pairing, unless you give paired runs the same sample name.")
  # conditions
  cond <- lapply(pm$sample, conditions_from_names, pattern = pattern)
  # a table whose names are bare condition names of the other table (e.g. SILAC labels L, H)
  for (k in names(cond)) {
    o <- setdiff(names(cond), k)
    if (all(is.na(cond[[k]])) && length(pm$sample[[k]]) && all(pm$sample[[k]] %in% cond[[o]])) cond[[k]] <- pm$sample[[k]]
  }
  cnd <- sort(unique(stats::na.omit(unlist(cond))))
  if (length(cnd)) {
    i <- which(!is.na(cond$ptm))[1]
    ex <- if (!is.na(i)) sprintf(", e.g. %s -> %s", pm$sample$ptm[i], cond$ptm[i]) else ""
    msg <- c(msg, sprintf("%d conditions read from the sample names (the name without its replicate number%s): %s.",
                          length(cnd), ex, short(cnd, 8)))
    few <- cnd[vapply(cnd, function(cc) min(sum(cond$protein %in% cc), sum(cond$ptm %in% cc)) < 2, logical(1))]
    if (length(few)) msg <- c(msg, sprintf("Check: fewer than 2 samples in %s.", short(few, 6)))
  } else msg <- c(msg, "Conditions could not be read from the sample names: type them in the tables or upload an annotation file.")
  ann <- lapply(c(protein = "protein", ptm = "ptm"), function(k) {
    d <- data.frame(run = runs[[k]], sample = "", condition = "", factor = "", subject = "", mixture = "", stringsAsFactors = FALSE)
    u <- which(use[[k]])
    d$sample[u] <- pm$sample[[k]]
    cc <- cond[[k]]; cc[is.na(cc)] <- ""
    d$condition[u] <- cc
    d[order(!use[[k]]), , drop = FALSE]                                       # excluded runs last
  })
  for (k in names(ann)) rownames(ann[[k]]) <- NULL
  list(annotation = ann, messages = msg)
}

#' Split counts into a high and a low group at a gap of more than `ratio`-fold
#' @return logical vector (TRUE = high group) named like `n`, or NULL when there is no such gap.
#' @keywords internal
gap_split <- function(n, ratio = 3) {
  v <- sort(log10(n + 1)); if (length(v) < 2) return(NULL)
  g <- diff(v); i <- which.max(g)
  if (g[i] < log10(ratio)) return(NULL)
  n > 10^v[i] - 1 + 1e-9
}

#' Tokens of run names, keeping the separators
#' @keywords internal
run_tokens <- function(s, pattern) {
  lapply(s, function(x) {
    tok <- regmatches(x, gregexpr(pattern, x, perl = TRUE), invert = TRUE)[[1]]
    sep <- c("", regmatches(x, gregexpr(pattern, x, perl = TRUE))[[1]])
    keep <- nzchar(tok); list(tok = tok[keep], sep = sep[seq_along(tok)][keep])
  })
}

#' Pair PTM and proteome runs by name
#'
#' Removes the parts of the names that only mark the fraction (present in every name of one
#' table and in less than half of the other's), pairs runs whose remaining names are equal, then
#' pairs the remaining runs with the most similar name when the best match is mutual, unique
#' and shares at least 60 % of the parts. Paired runs get the proteome run's remaining name.
#' @return list(sample = list(protein, ptm), removed = list(protein, ptm), n_fuzzy, fuzzy_example)
#' @keywords internal
pair_runs <- function(pr, pt, pattern = "[_ .-]") {
  out <- list(sample = list(protein = pr, ptm = pt), removed = list(protein = character(), ptm = character()),
              n_fuzzy = 0L, fuzzy_example = character())
  if (!length(pr) || !length(pt)) return(out)
  tk <- list(protein = run_tokens(pr, pattern), ptm = run_tokens(pt, pattern))
  low <- lapply(tk, function(z) lapply(z, function(t) tolower(t$tok)))
  in_all <- function(l) Reduce(intersect, lapply(l, unique))
  share <- function(t, l) mean(vapply(l, function(v) t %in% v, logical(1)))
  marker <- list(protein = Filter(function(t) share(t, low$ptm) < 0.5, in_all(low$protein)),
                 ptm = Filter(function(t) share(t, low$protein) < 0.5, in_all(low$ptm)))
  if (length(pr) < 2) marker$protein <- character(); if (length(pt) < 2) marker$ptm <- character()
  rebuild <- function(z, drop) vapply(z, function(t) {
    k <- !tolower(t$tok) %in% drop
    if (!any(k)) return(paste0(t$tok, collapse = "_"))
    tok <- t$tok[k]; sep <- t$sep[k]; sep[1] <- ""; paste0(sep, tok, collapse = "") }, "")
  nm <- list(protein = rebuild(tk$protein, marker$protein), ptm = rebuild(tk$ptm, marker$ptm))
  if (anyDuplicated(nm$protein) || anyDuplicated(nm$ptm)) {                # removal would merge runs: keep the full names
    nm <- list(protein = pr, ptm = pt); marker <- list(protein = character(), ptm = character())
  }
  out$removed <- lapply(names(marker), function(k) unique(unlist(lapply(tk[[k]], function(t) t$tok[tolower(t$tok) %in% marker[[k]]]))))
  names(out$removed) <- names(marker)
  key <- function(x) tolower(gsub(pattern, "_", x, perl = TRUE))
  smp_pt <- nm$ptm
  hit <- match(key(nm$ptm), key(nm$protein))
  smp_pt[!is.na(hit)] <- nm$protein[hit[!is.na(hit)]]
  # fuzzy pairing of the rest
  free_pt <- which(is.na(hit)); free_pr <- setdiff(seq_along(pr), hit)
  if (length(free_pt) && length(free_pr)) {
    tset <- function(x) unique(tolower(regmatches(x, gregexpr(pattern, x, perl = TRUE), invert = TRUE)[[1]]))
    A <- lapply(nm$ptm[free_pt], tset); B <- lapply(nm$protein[free_pr], tset)
    J <- outer(seq_along(A), seq_along(B), Vectorize(function(i, j) length(intersect(A[[i]], B[[j]])) / length(union(A[[i]], B[[j]]))))
    J <- matrix(J, length(A), length(B))
    for (i in seq_along(A)) {
      j <- which.max(J[i, ]); best <- J[i, j]
      if (best < 0.6 || sum(J[i, ] == best) > 1 || which.max(J[, j]) != i || sum(J[, j] == best) > 1) next
      smp_pt[free_pt[i]] <- nm$protein[free_pr[j]]
      out$n_fuzzy <- out$n_fuzzy + 1L
      if (!length(out$fuzzy_example)) out$fuzzy_example <- c(pt[free_pt[i]], pr[free_pr[j]])
    }
  }
  if (anyDuplicated(smp_pt)) smp_pt <- make.unique(smp_pt)
  out$sample <- list(protein = nm$protein, ptm = smp_pt)
  out
}

#' Read conditions from sample names
#'
#' The part that varies like a replicate number is removed (and a leading run number), the rest
#' is the condition; among the candidate positions, the one that gives every condition at least
#' two samples is chosen.
#' @return character vector of conditions (NA when no replicate structure is found).
#' @keywords internal
conditions_from_names <- function(s, pattern = "[_ .-]") {
  if (!length(s)) return(character())
  tok <- lapply(run_tokens(s, pattern), `[[`, "tok")
  repl <- "^(?:[Rr](?:ep)?|[Bb](?:io)?[Rr]?(?:ep)?|[Ss]|[Tt]ech|[Ii]nj)?[-_]?0*[0-9]{1,3}$"
  score <- function(cn) { if (anyNA(cn)) return(-1); t <- table(cn); if (length(t) < 2) return(-1); mean(cn %in% names(t)[t >= 2]) }
  cands <- list()
  k <- unique(lengths(tok))
  if (length(k) == 1 && k > 1) {
    pos <- which(vapply(seq_len(k), function(p) mean(grepl(repl, vapply(tok, `[`, "", p), perl = TRUE)) >= 0.8, logical(1)))
    lead <- if (k > 2 && all(grepl("^[0-9]+$", vapply(tok, `[`, "", 1)))) 1L else integer()
    const <- which(vapply(seq_len(k), function(p) length(unique(vapply(tok, `[`, "", p))) == 1, logical(1)))   # e.g. "run", "QE"
    for (p in pos) for (drop in unique(list(p, unique(c(lead, p))))) {
      if (length(drop) >= k) next
      drop2 <- if (length(setdiff(seq_len(k), union(drop, const)))) union(drop, const) else drop
      cn <- vapply(tok, function(t) paste(t[-drop2], collapse = "_"), "")
      cands[[length(cands) + 1]] <- list(cn = cn, s = score(cn), p = max(drop))
    }
  }
  # names of different lengths: the condition is what precedes the last replicate-like part
  cn <- vapply(tok, function(t) {
    if (length(t) > 2 && grepl("^[0-9]+$", t[1])) t <- t[-1]
    r <- which(grepl(repl, t, perl = TRUE)); r <- r[r > 1]
    end <- if (length(r)) max(r) else length(t)
    if (end < 2) NA_character_ else paste(t[seq_len(end - 1)], collapse = "_") }, "")
  cands[[length(cands) + 1]] <- list(cn = cn, s = score(cn), p = -1)
  sc <- vapply(cands, `[[`, 0, "s")
  if (max(sc) < 0.5) return(rep(NA_character_, length(s)))
  best <- which(sc == max(sc)); best <- best[which.max(vapply(cands[best], `[[`, 0, "p"))]
  cands[[best]]$cn
}

#' Rename the runs (columns) of an assay
#' @param a a `ptm_assay`.
#' @param new new names, in the order of the columns (or named by the old names).
#' @return the assay with renamed columns.
#' @export
set_run_names <- function(a, new) {
  old <- colnames(a$intensity)
  if (!is.null(names(new))) new <- unname(new[old])
  if (length(new) != length(old) || anyNA(new) || anyDuplicated(new)) stop("New run names must be unique and one per run")
  colnames(a$intensity) <- new
  if (!is.null(a$raw_detected)) colnames(a$raw_detected) <- new
  if (!is.null(a$n_features)) colnames(a$n_features) <- new
  a
}
