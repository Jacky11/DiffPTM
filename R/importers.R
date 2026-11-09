#' Column specifications for supported quantification software
#'
#' Each spec describes where to find identifiers, quantitative columns and flags in the
#' protein and PTM tables of a given software. Specs are plain lists so they can be
#' overridden from the app ("Advanced: column mapping") without touching the code.
#'
#' Elements:
#' * `protein`: `id` (column(s) tried in order), `gene`, `quant` (one or several regex on
#'   column names tried in order, the first with matches wins; capture group 1 = sample name), `flags` (columns whose non-empty value flags a row
#'   for removal, e.g. reverse hits), `sep`.
#' * `ptm`: `id` (optional; built if absent), `protein_ids`, `gene`, `stripped_sequence`,
#'   `modified_sequence`, `residue`, `position`, `localization_prob` (regex, several columns
#'   allowed), `ptm_type`, `multiplicity` (column giving the number of modifications on the
#'   peptide), `multiplicity_suffix` (regex removed from quant column names, e.g. `___[123]$`),
#'   `quant`, `flags`, `sep`. Column names may be given as a vector of alternatives tried in
#'   order.
#'
#' @param software one of `names(column_specs())`, or `NULL` for the whole registry.
#' @export
column_specs <- function(software = NULL) {
  specs <- list(
    MaxQuant = list(
      protein = list(
        id = c("Majority protein IDs", "Protein IDs"),
        gene = "Gene names",
        quant = c("^LFQ intensity (.+)$", "^Reporter intensity corrected (.+)$", "^Intensity (.+)$"),
        flags = c("Reverse", "Potential contaminant", "Only identified by site"),
        sep = ";"
      ),
      ptm = list(
        protein_ids = "Proteins",
        gene = "Gene names",
        stripped_sequence = NULL,          # Phospho(STY)Sites.txt has no stripped sequence
        modified_sequence = "Sequence window",
        residue = "Amino acid",
        position = "Position",
        positions_within_proteins = "Positions within proteins",
        localization_prob = "^Localization prob$",
        ptm_type = NULL,
        multiplicity_suffix = "___([123])$",
        quant = c("^Reporter intensity corrected (.+?)(?:___[123])?$", "^Intensity (.+?)(?:___[123])?$"),
        flags = c("Reverse", "Potential contaminant"),
        sep = ";"
      )
    ),
    Spectronaut = list(                       # pivot reports (semicolon- or tab-separated, decimal comma tolerated)
      protein = list(
        id = c("PG.ProteinGroups", "PG.ProteinAccessions"),
        gene = "PG.Genes",
        quant = "^(?:\\[\\d+\\] )?(.+?)\\.(?:PG\\.Quantity|PG\\.MS2Quantity)$",
        flags = character(0),
        sep = ";"
      ),
      ptm = list(                             # PTM site report: one row per site x protein accession
        id = "PTM.CollapseKey",
        protein_ids = "PTM.ProteinId",
        gene = "PG.Genes",
        stripped_sequence = "PEP.StrippedSequence",
        modified_sequence = c("PTM.FlankingRegion", "PTM.CollapseKey"),
        residue = "PTM.SiteAA",
        position = "PTM.SiteLocation",
        localization_prob = "PTM\\.SiteProbability$",
        ptm_type = "PTM.ModificationTitle",
        multiplicity = "PTM.Multiplicity",
        multiplicity_suffix = NULL,
        quant = "^(?:\\[\\d+\\] )?(.+?)\\.PTM\\.Quantity$",
        flags = character(0),
        sep = ";"
      )
    ),
    `Proteome Discoverer` = list(             # Proteins and PeptideGroups exports (tab-separated text)
      protein = list(
        id = c("Accession", "Master Protein Accessions"),
        gene = c("Gene Symbol", "Gene ID"),
        # per-file abundances: "Abundance: F1: Sample, control" (LFQ) or "Abundance: F1: 126, control" (TMT);
        # normalised values preferred; grouped / scaled / ratio columns are not per-sample and are ignored
        quant = c("^Abundances \\(Normalized\\): (F\\d+: [^,]+)(?:,.*)?$", "^Abundance: (F\\d+: [^,]+)(?:,.*)?$",
                  "^Abundances \\(Scaled\\): (F\\d+: [^,]+)(?:,.*)?$", "^Abundance (F\\d+)(?: .*)?$"),
        flags = c("Contaminant"),
        sep = "; "
      ),
      ptm = list(                             # one row per peptide group (peptidoform); sites parsed from
        protein_ids = "Master Protein Accessions",   # "Modifications in Master Proteins"
        gene = c("Gene Symbol", "Master Protein Gene Symbol"),
        stripped_sequence = c("Sequence", "Annotated Sequence"),
        modified_sequence = "Modifications in Master Proteins",
        residue = NULL,
        position = NULL,
        localization_prob = NULL,
        ptm_type = NULL,
        multiplicity_suffix = NULL,
        quant = c("^Abundances \\(Normalized\\): (F\\d+: [^,]+)(?:,.*)?$", "^Abundance: (F\\d+: [^,]+)(?:,.*)?$",
                  "^Abundances \\(Scaled\\): (F\\d+: [^,]+)(?:,.*)?$", "^Abundance (F\\d+)(?: .*)?$"),
        flags = c("Contaminant"),
        sep = "; "
      )
    ),
    FragPipe = list(
      protein = list(                         # combined_protein.tsv (label-free) or abundance_protein_*.tsv (TMT-Integrator)
        id = c("Protein ID", "ProteinID", "Index", "Protein"),
        gene = "Gene",
        quant = c("^(.+?) MaxLFQ Intensity$", "^(.+?) Intensity$",
                  "^(?!Index$|Gene$|ProteinID$|Protein$|Peptide$|SequenceWindow$|Start$|End$|MaxPepProb$|ReferenceIntensity$|NumberPSM$|Entry Name$|Description$|Modifications$)(.+)$"),
        flags = character(0),
        sep = ","
      ),
      ptm = list(                             # abundance_single-site_*.tsv (TMT-Integrator) or combined_modified_peptide.tsv (label-free)
        id = "Index",
        protein_ids = c("ProteinID", "Protein ID"),
        gene = "Gene",
        stripped_sequence = c("Peptide", "Peptide Sequence"),
        modified_sequence = c("Modified Sequence", "Index"),
        residue = NULL,
        position = NULL,                      # parsed from Index "P12345_S123" or from Assigned Modifications by import_ptm()
        localization_prob = c("^MaxPepProb$", "^Probability$"),
        ptm_type = NULL,
        multiplicity_suffix = NULL,
        quant = c("^(.+?) MaxLFQ Intensity$", "^(.+?) Intensity$",
                  "^(?!Index$|Gene$|ProteinID$|Protein$|Peptide$|SequenceWindow$|Start$|End$|MaxPepProb$|ReferenceIntensity$|NumberPSM$|Entry Name$|Description$|Modifications$)(.+)$"),
        flags = character(0),
        sep = ","
      )
    ),
    `DIA-NN` = list(                          # sample columns are named by the raw file path
      protein = list(                         # report.pg_matrix.tsv
        id = c("Protein.Group", "Protein.Ids"),
        gene = "Genes",
        quant = "^(?!Protein\\.Group$|Protein\\.Ids$|Protein\\.Names$|Genes$|First\\.Protein\\.Description$|N\\.Sequences$|N\\.Proteotypic\\.Sequences$)(.+)$",
        flags = character(0),
        sep = ";"
      ),
      ptm = list(                             # report.pr_matrix.tsv (precursors) or a site matrix (report.phosphosites_*.tsv)
        protein_ids = c("Protein.Group", "Protein.Ids"),
        gene = "Genes",
        stripped_sequence = "Stripped.Sequence",
        modified_sequence = c("Modified.Sequence", "Site.ID", "Precursor.Id"),
        residue = NULL,
        position = c("Site.Position", "Position"),
        localization_prob = "^(?:PTM\\.Site\\.Confidence|Site\\.Confidence|Site\\.Localization\\.Confidence|PTM\\.Q\\.Value)$",
        ptm_type = NULL,
        multiplicity_suffix = NULL,
        quant = "^(?!Protein\\.Group$|Protein\\.Ids$|Protein\\.Names$|Genes$|First\\.Protein\\.Description$|Proteotypic$|Stripped\\.Sequence$|Modified\\.Sequence$|Precursor\\.Charge$|Precursor\\.Id$|Site\\.ID$|Site\\.Position$|Position$|Site\\.Confidence$|PTM\\.Site\\.Confidence$|Site\\.Localization\\.Confidence$|PTM\\.Q\\.Value$|N\\.Sequences$|N\\.Proteotypic\\.Sequences$)(.+)$",
        flags = character(0),
        sep = ";"
      )
    ),
    Generic = list(                           # wide table: id columns + sample columns
      protein = list(id = c("protein_id", "feature_id", "protein_ids"), gene = "gene", quant = "^intensity[_. ](.+)$",
                     flags = character(0), sep = ";"),
      ptm = list(id = "feature_id", protein_ids = "protein_ids", gene = "gene",
                 stripped_sequence = "stripped_sequence", modified_sequence = "modified_sequence",
                 residue = "residue", position = "position", localization_prob = "^localization_prob$",
                 ptm_type = "ptm_type", multiplicity_suffix = NULL,
                 quant = "^intensity[_. ](.+)$", flags = character(0), sep = ";")
    )
  )
  if (is.null(software)) return(specs)
  if (!software %in% names(specs)) stop("Unknown software: ", software,
                                        ". Known: ", paste(names(specs), collapse = ", "))
  specs[[software]]
}

#' Read a tsv / csv / txt / xlsx table without mangling column names
#' @param path file path.
#' @param sheet sheet index or name for Excel files.
#' @param dec decimal mark, `"."` or `","`; guessed from the file when `NULL`.
#' @export
import_table <- function(path, sheet = 1, dec = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Package readxl is needed for Excel input")
    return(as.data.frame(readxl::read_excel(path, sheet = sheet), stringsAsFactors = FALSE))
  }
  na <- c("", "NA", "NaN", "Filtered", "#N/A", "n. def.", "n.def.")
  if (is.null(dec)) dec <- guess_decimal_mark(path, na)
  data.table::fread(path, data.table = FALSE, check.names = FALSE, na.strings = na, dec = dec, encoding = "UTF-8", integer64 = "double")
}

#' Guess whether a text table uses the decimal comma
#'
#' Spectronaut and Proteome Discoverer exports made on a computer with a European locale
#' write `60441,39`; the first 2,000 rows are scanned and `","` is returned when at least one
#' column is made of such numbers.
#' @param path file path.
#' @param na strings read as missing values.
#' @return `"."` or `","`.
#' @keywords internal
guess_decimal_mark <- function(path, na = c("", "NA", "NaN", "Filtered")) {
  d <- tryCatch(data.table::fread(path, data.table = FALSE, check.names = FALSE, na.strings = na, nrows = 2000,
                                  colClasses = "character", encoding = "UTF-8"), error = function(e) NULL)
  if (is.null(d) || !ncol(d)) return(".")
  comma <- "^-?[0-9]+,[0-9]+(?:[eE][-+]?[0-9]+)?$"
  for (v in d) {
    v <- v[!is.na(v)]
    if (length(v) >= 2 && mean(grepl(comma, v, perl = TRUE)) > 0.5 && !any(grepl(",.*,", v))) return(",")
  }
  "."
}

#' Warn when the scale of the values does not match `log2_transform`
#'
#' Raw intensities are large (1e3-1e10); values that never exceed 64 are almost surely on
#' the log2 scale already (TMT-Integrator, MSstats summaries), and values above 1e4 are not.
#' @keywords internal
check_scale <- function(m, log2_transform) {
  v <- m[is.finite(m)]
  if (!length(v)) return(log2_transform)
  if (log2_transform && max(v) < 64 && min(v) > -64)
    warning("Values never exceed 64: they look log2-transformed already; the log2 transformation was NOT applied (untick 'Intensities are raw' / use log2_transform = FALSE to silence this).", call. = FALSE)
  if (!log2_transform && max(v) > 1e4)
    warning("Values exceed 10,000: they look like raw intensities, but log2_transform = FALSE.", call. = FALSE)
  if (log2_transform && max(v) < 64 && min(v) > -64) FALSE else log2_transform
}

first_present <- function(cols, d) {
  cols <- cols[!vapply(cols, is.null, TRUE)]
  hit <- cols[cols %in% names(d)]
  if (length(hit)) hit[1] else NULL
}

#' Extract the quantitative matrix from a wide table using a regex with one capture group
#' @keywords internal
quant_matrix <- function(d, pattern, suffix = NULL) {
  cols <- character(0)
  for (p in pattern) { cols <- grep(p, names(d), value = TRUE, perl = TRUE); if (length(cols)) { pattern <- p; break } }
  if (!length(cols)) stop("No quantitative columns matched the pattern(s): ", paste(pattern, collapse = " | "))
  samples <- sub(pattern, "\\1", cols, perl = TRUE)
  m <- as.matrix(d[, cols, drop = FALSE])
  m <- suppressWarnings(apply(m, 2, as.numeric))
  if (is.null(dim(m))) m <- matrix(m, ncol = length(cols))
  colnames(m) <- if (is.null(suffix)) samples else cols
  attr(m, "pattern") <- pattern
  m
}

#' Import a total-proteome table
#'
#' @param path file path (or a data.frame already read).
#' @param software name in [column_specs()].
#' @param spec optional spec overriding `column_specs(software)$protein`.
#' @param quant_pattern optional regex overriding `spec$quant`.
#' @param log2_transform transform raw intensities (default `TRUE`; set `FALSE` for
#'   tables already on log scale).
#' @return a `ptm_assay`.
#' @export
import_protein <- function(path, software = "MaxQuant", spec = NULL, quant_pattern = NULL,
                           log2_transform = TRUE) {
  d <- if (is.data.frame(path)) path else import_table(path)
  spec <- modifyList(column_specs(software)$protein, spec %||% list())
  if (!is.null(quant_pattern)) spec$quant <- quant_pattern
  idcol <- first_present(spec$id, d)
  if (is.null(idcol)) stop("None of the identifier columns found: ", paste(spec$id, collapse = ", "))
  m <- quant_matrix(d, spec$quant)
  log2_transform <- check_scale(m, log2_transform)
  meta <- data.frame(feature_id = make.unique(as.character(d[[idcol]])),
                     protein_ids = as.character(d[[idcol]]),
                     gene = if (!is.null(first_present(spec$gene, d))) as.character(d[[first_present(spec$gene, d)]]) else NA_character_,
                     stringsAsFactors = FALSE)
  for (f in spec$flags) if (f %in% names(d)) meta[[make.names(f)]] <- !is.na(d[[f]]) & nzchar(as.character(d[[f]])) & as.character(d[[f]]) != "FALSE"
  colnames(m) <- clean_sample_names(colnames(m), software)
  g <- ifelse(is.na(meta$gene) | !nzchar(meta$gene), meta$protein_ids, meta$gene)
  meta$label <- vapply(strsplit(g, spec$sep, fixed = TRUE), function(t) if (length(t)) trimws(t[1]) else "", "")
  new_assay(m, meta, log2_transform = log2_transform, sep = spec$sep)
}

#' Import a PTM table (sites or modified peptides)
#'
#' @inheritParams import_protein
#' @param multiplicity how to handle MaxQuant `___1/2/3` columns: `"expand"` creates one
#'   feature per multiplicity (as the original app's checkbox), `"first"` keeps `___1` only,
#'   `"sum"` sums the three on the linear scale.
#' @param use_max_localization when several localization-probability columns match (e.g.
#'   Spectronaut per-sample columns), keep the maximum per feature.
#' @export
import_ptm <- function(path, software = "MaxQuant", spec = NULL, quant_pattern = NULL,
                       log2_transform = TRUE, multiplicity = c("expand", "first", "sum"),
                       use_max_localization = TRUE) {
  multiplicity <- match.arg(multiplicity)
  d <- if (is.data.frame(path)) path else import_table(path)
  spec <- modifyList(column_specs(software)$ptm, spec %||% list())
  if (!is.null(quant_pattern)) spec$quant <- quant_pattern

  getcol <- function(name) {
    cn <- first_present(spec[[name]], d)
    if (is.null(cn)) return(rep(NA_character_, nrow(d)))
    as.character(d[[cn]])
  }
  loc <- rep(NA_real_, nrow(d))
  if (!is.null(spec$localization_prob)) {
    lc <- character(0)
    for (lp in spec$localization_prob) { lc <- grep(lp, names(d), value = TRUE, perl = TRUE); if (length(lc)) break }
    if (length(lc)) {
      lm <- suppressWarnings(apply(as.matrix(d[, lc, drop = FALSE]), 2, as.numeric))
      if (is.null(dim(lm))) lm <- matrix(lm, ncol = length(lc))
      loc <- if (use_max_localization) apply(lm, 1, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)) else lm[, 1]
    }
  }
  meta <- data.frame(
    protein_ids = getcol("protein_ids"),
    gene = getcol("gene"),
    stripped_sequence = getcol("stripped_sequence"),
    modified_sequence = getcol("modified_sequence"),
    residue = getcol("residue"),
    position = getcol("position"),
    positions_within_proteins = getcol("positions_within_proteins"),
    ptm_type = getcol("ptm_type"),
    localization_prob = loc,
    stringsAsFactors = FALSE
  )
  # FragPipe single-site index "P12345_S123" -> residue / position; combined_modified_peptide.tsv
  # ("Assigned Modifications" = "5S(79.9663), 9T(79.9663)", positions within the peptide + "Start")
  if (software == "FragPipe" && all(is.na(meta$position))) {
    if ("Index" %in% names(d) && all(grepl("_[A-Z]\\d+$", d$Index))) {
      site <- sub("^.*_([A-Z])(\\d+)$", "\\1\\2", as.character(d$Index))
      meta$residue  <- substr(site, 1, 1)
      meta$position <- sub("^[A-Z]", "", site)
    } else if (all(c("Assigned Modifications", "Start") %in% names(d))) {
      am <- regmatches(as.character(d[["Assigned Modifications"]]), gregexpr("(\\d+)([A-Z])\\(([0-9.]+)\\)", as.character(d[["Assigned Modifications"]]), perl = TRUE))
      start <- as.integer(d$Start)
      sites <- mapply(function(v, st) {
        if (!length(v)) return(c(NA, NA, NA))
        v <- v[!grepl("\\((?:57\\.0215|15\\.9949|42\\.0106|229\\.1629|304\\.2071)\\)", v, perl = TRUE)]   # Carbamidomethyl, Oxidation, Acetyl, TMT
        if (!length(v)) return(c(NA, NA, NA))
        c(paste(sub("^\\d+([A-Z]).*$", "\\1", v), collapse = ";"),
          paste(as.integer(sub("^(\\d+).*$", "\\1", v)) + st - 1L, collapse = ";"),
          paste(unique(sub("^.*\\(([0-9.]+)\\)$", "\\1", v)), collapse = ";"))
      }, am, start)
      meta$residue  <- sites[1, ]; meta$position <- sites[2, ]; meta$ptm_type <- sites[3, ]
    }
  }
  if (software == "Proteome Discoverer") meta <- parse_pd_modifications(meta, d)
  for (f in spec$flags) if (f %in% names(d)) meta[[make.names(f)]] <- !is.na(d[[f]]) & nzchar(as.character(d[[f]])) & as.character(d[[f]]) != "FALSE"

  m <- quant_matrix(d, spec$quant, suffix = spec$multiplicity_suffix)
  log2_transform <- check_scale(m, log2_transform)
  used_pattern <- attr(m, "pattern")
  idcol <- first_present(spec$id, d)
  first_acc <- vapply(split_ids(meta$protein_ids, spec$sep, strip_isoforms = FALSE), function(v) if (length(v)) v[1] else NA_character_, "")
  has_site <- !is.na(meta$residue) & !is.na(meta$position) & nzchar(meta$position) & !is.na(first_acc)
  base_id <- if (!is.null(idcol)) as.character(d[[idcol]]) else
    ifelse(has_site, paste0(first_acc, "_", meta$residue, meta$position),              # P12345_S123 (FragPipe-like)
           paste(meta$protein_ids, meta$residue, meta$position, meta$modified_sequence, sep = "|"))
  meta$feature_id <- make.unique(base_id)

  if (!is.null(spec$multiplicity_suffix) && any(grepl(spec$multiplicity_suffix, colnames(m), perl = TRUE))) {
    out <- expand_multiplicity(m, meta, spec$multiplicity_suffix, used_pattern, mode = multiplicity)
    m <- out$intensity; meta <- out$meta
  } else {
    colnames(m) <- sub(used_pattern, "\\1", colnames(m), perl = TRUE)
    meta$multiplicity <- getcol("multiplicity")
  }
  colnames(m) <- clean_sample_names(colnames(m), software)
  meta$label <- site_label(meta, sep = spec$sep)
  new_assay(m, meta, log2_transform = log2_transform, sep = spec$sep)
}

#' Short display label of a PTM feature: gene (or first accession), residue and position
#'
#' `ESYT2 S663` for a phosphosite, `ESYT2 S663 (x2)` when the feature is the doubly modified
#' form, `ESYT2 AAAsPAK` when no position is known.
#' @param meta feature metadata.
#' @param sep separator of protein identifiers.
#' @keywords internal
site_label <- function(meta, sep = ";") {
  first <- function(v) vapply(strsplit(ifelse(is.na(v), "", v), sep, fixed = TRUE),
                              function(t) if (length(t)) sub("^(?:sp|tr)\\|([^|]+)\\|.*$", "\\1", trimws(t[1]), perl = TRUE) else "", "")
  who <- first(meta$gene); acc <- first(meta$protein_ids)
  who[!nzchar(who)] <- acc[!nzchar(who)]
  site <- paste0(ifelse(is.na(meta$residue), "", meta$residue), ifelse(is.na(meta$position), "", meta$position))
  site[!nzchar(site)] <- ifelse(is.na(meta$modified_sequence), "", meta$modified_sequence)[!nzchar(site)]
  lab <- trimws(paste(who, site))
  mult <- suppressWarnings(as.integer(meta$multiplicity))
  lab <- ifelse(!is.na(mult) & mult > 1, sprintf("%s (x%d)", lab, mult), lab)
  lab[!nzchar(trimws(lab))] <- meta$feature_id[!nzchar(trimws(lab))]
  lab
}

#' Tidy sample names produced by a given software
#'
#' DIA-NN names columns by the full path of the raw file; Proteome Discoverer by
#' `F1: Sample` (label-free) or `F1: 126` (TMT); the `.raw` / `.d` extension and the
#' `Sample` channel are removed.
#' @keywords internal
clean_sample_names <- function(x, software) {
  if (software == "DIA-NN") x <- sub("\\.(?:raw|d|mzML|wiff|dia)$", "", basename(gsub("\\\\", "/", x)), ignore.case = TRUE, perl = TRUE)
  if (software == "Proteome Discoverer") x <- sub(": Sample$", "", x)
  make.unique(x)
}

#' Parse Proteome Discoverer "Modifications in Master Proteins"
#'
#' `P12345 1xPhospho [S123(100); T127(98.5)]` gives residue `S;T`, position `123;127`,
#' localization probability 0.985 (the lowest site probability, on a 0-1 scale) and
#' `ptm_type` `Phospho`. Carbamidomethyl, Oxidation and TMT labels are ignored.
#' @keywords internal
parse_pd_modifications <- function(meta, d) {
  src <- meta$modified_sequence
  if (all(is.na(src))) return(meta)
  skip <- "Carbamidomethyl|Oxidation|TMT|Acetyl \\[N-Term\\]|Met-loss|iTRAQ|Label"
  parsed <- lapply(src, function(s) {
    if (is.na(s)) return(list(res = NA, pos = NA, prob = NA, type = NA))
    mods <- regmatches(s, gregexpr("\\d+x[A-Za-z0-9 ()-]+? \\[[^]]*\\]", s, perl = TRUE))[[1]]
    mods <- mods[!grepl(skip, mods, perl = TRUE)]
    if (!length(mods)) return(list(res = NA, pos = NA, prob = NA, type = NA))
    type <- unique(sub("^\\d+x([A-Za-z0-9 ()-]+?) \\[.*$", "\\1", mods, perl = TRUE))
    sites <- unlist(regmatches(mods, gregexpr("[A-Z]\\d+(?:\\([0-9.]+\\))?", mods, perl = TRUE)))
    sites <- sites[!grepl("^N-Term|^C-Term", sites)]
    res <- substr(sites, 1, 1); pos <- sub("^[A-Z](\\d+).*$", "\\1", sites)
    pr <- suppressWarnings(as.numeric(sub("^[A-Z]\\d+\\(([0-9.]+)\\)$", "\\1", sites)))
    list(res = paste(res, collapse = ";"), pos = paste(pos, collapse = ";"),
         prob = if (all(is.na(pr))) NA_real_ else min(pr, na.rm = TRUE) / 100, type = paste(type, collapse = ";"))
  })
  meta$residue  <- vapply(parsed, function(p) as.character(p$res), "")
  meta$position <- vapply(parsed, function(p) as.character(p$pos), "")
  pr <- vapply(parsed, function(p) as.numeric(p$prob), 0)
  meta$localization_prob <- ifelse(is.na(meta$localization_prob), pr, meta$localization_prob)
  tp <- vapply(parsed, function(p) as.character(p$type), "")
  meta$ptm_type <- ifelse(is.na(tp), meta$ptm_type, tp)
  if ("Positions in Master Proteins" %in% names(d) && all(is.na(meta$position) | !nzchar(meta$position))) {
    meta$position <- sub("^.*\\[(\\d+)-\\d+\\].*$", "\\1", as.character(d[["Positions in Master Proteins"]]))
  }
  meta$stripped_sequence <- toupper(gsub("^\\[[^]]*\\]\\.|\\.\\[[^]]*\\]$", "", meta$stripped_sequence))
  meta
}

#' Handle MaxQuant multiplicity columns (`Intensity sample___1/2/3`)
#'
#' @param m raw (not log) matrix whose column names still carry the suffix.
#' @param meta feature metadata.
#' @param suffix regex with the multiplicity number as capture group 1.
#' @param quant_pattern regex giving the sample name as capture group 1.
#' @param mode `"expand"`, `"first"` or `"sum"`.
#' @return list(intensity, meta) with a `multiplicity` column in `meta`.
#' @export
expand_multiplicity <- function(m, meta, suffix = "___([123])$", quant_pattern, mode = "expand") {
  m <- m[, grepl(suffix, colnames(m), perl = TRUE), drop = FALSE]   # drop "Intensity X" totals
  mult   <- sub(paste0("^.*", suffix), "\\1", colnames(m), perl = TRUE)
  sample <- sub(quant_pattern, "\\1", colnames(m), perl = TRUE)
  levels_m <- sort(unique(mult))
  by_m <- lapply(levels_m, function(k) {
    sub_m <- m[, mult == k, drop = FALSE]; colnames(sub_m) <- sample[mult == k]; sub_m
  })
  names(by_m) <- levels_m
  if (mode == "first") {
    meta$multiplicity <- levels_m[1]
    return(list(intensity = by_m[[1]], meta = meta))
  }
  if (mode == "sum") {
    s <- Reduce(function(a, b) { a[is.na(a)] <- 0; b[is.na(b)] <- 0; a + b }, by_m)
    s[s == 0] <- NA
    meta$multiplicity <- "sum"
    return(list(intensity = s, meta = meta))
  }
  # expand: one feature per multiplicity
  int <- do.call(rbind, by_m)
  meta_x <- do.call(rbind, lapply(levels_m, function(k) { mm <- meta; mm$multiplicity <- k; mm }))
  meta_x$feature_id <- paste0(meta_x$feature_id, "___", meta_x$multiplicity)
  keep <- rowSums(!is.na(int) & int != 0) > 0        # raw scale: 0 means not observed
  list(intensity = int[keep, , drop = FALSE], meta = meta_x[keep, , drop = FALSE])
}

`%||%` <- function(a, b) if (is.null(a)) b else a
