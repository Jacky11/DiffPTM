#################################################
# Reading of the input file (MaxQuant proteinGroups.txt or Spectronaut report)
#################################################

# Convert a Spectronaut report to the MaxQuant column naming used by the pipeline
spectronaut_convert <- function(data) {

  # Process header names
  header = names(data)
  header[grep("PG.ProteinGroups", header)] <- "Protein.IDs"
  header[grep("PG.Genes", header)] <- "Majority.protein.IDs"
  if (length(grep("PG.FastaHeaders", header)) > 0) {
    header[grep("PG.FastaHeaders", header)] <- "Fasta.headers"
  } else {
    header[grep("PG.ProteinAccessions", header)] <- "Fasta.headers"
  }
  if (length(grep("PG.ProteinDescriptions", header)) > 0) {
    header[grep("PG.ProteinDescriptions", header)] <- "Protein.Description"
  }

  header = gsub("^X\\.\\d+\\.\\.", "", header)

  # Number of stripped sequences identified -> peptides / unique peptides
  id_subh = grep(".PG.NrOfStrippedSequencesIdentified", header)
  dUnique = data[, id_subh]

  prefix_commun = find_common_prefix(header[id_subh])
  diff = sub(prefix_commun, "", header[id_subh])
  diff = reverse_string(diff)
  suffix_commun = find_common_prefix(diff)
  diff = sub(suffix_commun, "", diff)
  diff = reverse_string(diff)

  cdu = diff
  for (i in 1:length(cdu)) {
    nn = unlist(strsplit(x = cdu[i], split = "[.]"))
    if (length(nn) > 3) {
      cdu[i] = paste0(nn[1:(length(nn) - 3)], collapse = ".")
    }
  }
  diff = cdu

  num_diff = suppressWarnings(as.numeric(diff))
  isNumeric = all(!is.na(num_diff))

  if (isNumeric) {
    header[id_subh] <- paste("Peptides.Smple", diff, sep = ".")
    names(dUnique) <- paste("Unique.peptides.Smple", diff, sep = ".")
  } else {
    header[id_subh] <- paste("Peptides", diff, sep = ".")
    names(dUnique) <- paste("Unique.peptides", diff, sep = ".")
  }

  # Quantity columns
  id_subh <- grep("PG.Quantity", header)
  prefix_commun = find_common_prefix(header[id_subh])
  header_differences = sub(prefix_commun, "", header[id_subh])
  header_differences = reverse_string(header_differences)
  suffix_commun = find_common_prefix(header_differences)
  header_differences = sub(suffix_commun, "", header_differences)
  header_differences = reverse_string(header_differences)

  if (isNumeric) {
    header[id_subh] <- paste("Intensity.Smple", header_differences, sep = ".")
  } else {
    header[id_subh] <- paste("Intensity", header_differences, sep = ".")
  }
  data[, id_subh] <- apply(data[, id_subh], 2, convert_number)

  # iBAQ columns, if present
  indexGrep = grep("PG.IBAQ", header)
  if (length(indexGrep) > 0) {
    id_subh <- grep("PG.IBAQ", header)
    prefix_commun = find_common_prefix(header[id_subh])
    header_differences = sub(prefix_commun, "", header[id_subh])
    header_differences = reverse_string(header_differences)
    suffix_commun = find_common_prefix(header_differences)
    header_differences = sub(suffix_commun, "", header_differences)
    header_differences = reverse_string(header_differences)

    if (isNumeric) {
      header[id_subh] <- paste("iBAQ.Smple", header_differences, sep = ".")
    } else {
      header[id_subh] <- paste("iBAQ", header_differences, sep = ".")
    }

    # semicolon-delimited values: keep the maximum
    for (col in id_subh) {
      data[, col] <- sapply(data[, col], function(cell) {
        if (contains_semicolon(cell)) find_max(cell) else convert_number(cell)
      })
    }
    data[, id_subh] <- as.data.frame(data[, id_subh])
    iBAQ <- apply(data[, id_subh], 1, function(x) sum(x, na.rm = TRUE))
  }

  # Placeholder columns for MaxQuant compatibility
  OnlyidBySite <- rep(" ", nrow(data))
  Reverse <- rep(" ", nrow(data))
  PotentialContaminants <- rep(" ", nrow(data))

  colnames(data) <- header
  if (length(indexGrep) > 0) {
    final <- cbind(data, iBAQ, dUnique, OnlyidBySite, Reverse, PotentialContaminants)
  } else {
    final <- cbind(data, dUnique, OnlyidBySite, Reverse, PotentialContaminants)
  }
  colnames(final)[grep("OnlyidBySite", colnames(final))] <- "Only.identified.by.site"
  colnames(final)[grep("PotentialContaminant", colnames(final))] <- "Potential.contaminant"
  colnames(final)[grep("Reverse", colnames(final))] <- "Reverse"

  final[is.na(final)] <- 0
  return(final)
}


# Read the input file and extract the matrices used by the pipeline.
#   file        : path to proteinGroups.txt (MaxQuant) or the Spectronaut report
#   spectronaut : TRUE for a Spectronaut report
#   intensity   : "Raw" (Intensity columns) or "LFQ" (LFQ intensity columns)
#   name, size  : file name / size reported on the first slide (Shiny passes the
#                 name of the uploaded file, which differs from its temporary path)
read_input <- function(file, spectronaut = FALSE, intensity = c("Raw", "LFQ"),
                       name = basename(file), size = file.size(file)) {
  intensity <- match.arg(intensity)

  dat <- NULL
  if (spectronaut) {
    # separator detected on the header line (tab, semicolon or comma)
    first <- readLines(file, n = 1, warn = FALSE)
    counts <- c(tab = lengths(regmatches(first, gregexpr("\t", first))),
                semicolon = lengths(regmatches(first, gregexpr(";", first))),
                comma = lengths(regmatches(first, gregexpr(",", first))))
    dat <- switch(names(which.max(counts)),
                  tab = read.csv(file, header = TRUE, sep = "\t", quote = ""),
                  semicolon = read.csv(file, header = TRUE, sep = ";"),
                  comma = read.csv(file, header = TRUE, sep = ","))
    dat = spectronaut_convert(dat)
  } else {
    dat = read.csv(file, header = TRUE, sep = "\t", quote = "")
  }

  # Intensities
  if (intensity == "LFQ") {
    rdat = as.matrix(dat[, grep("LFQ", names(dat), value = TRUE)])
  }
  if (intensity == "Raw") {
    rdat = as.matrix(dat[, grep("Intensity.", names(dat), value = TRUE)])
  }

  if ("Protein.Ids" %in% names(dat)) {
    colnames(dat)[grep("Protein.Ids", names(dat))] = "Protein.IDs"
  }

  # Peptide counts
  rdat_id = data.frame(as.matrix(dat[, grep("Peptides.", names(dat), value = TRUE)]))
  rdat_id2 = data.frame(as.matrix(dat[, grep("Unique.peptides.", names(dat), value = TRUE)]))

  # Decimal commas
  if (is.character(rdat[, 1])) { rdat = gsub(",", ".", rdat) }

  # iBAQ columns containing semicolon-separated values: keep the maximum
  colonnes_iBAQ <- grep("iBAQ", colnames(dat))
  if (length(colonnes_iBAQ) > 0) {
    contient_point_virgule <- any(sapply(dat[, colonnes_iBAQ], function(x) any(grepl(";", x))))
    if (contient_point_virgule) {
      dat[, colonnes_iBAQ] <- apply(dat[, colonnes_iBAQ], 2, function(x) {
        sapply(x, function(val) {
          if (grepl(";", val)) {
            val = gsub('"', '', val)
            valeur_num <- unlist(lapply(strsplit(val, ";")[[1]], convert_number))
            valeur_num = max(valeur_num, na.rm = TRUE)
          } else {
            if (is.na(val)) return(NA)
            valeur_num = convert_number(val)
          }
          return(valeur_num)
        })
      })
    }
  }

  rdat = apply(rdat, 2, as.numeric)
  rdat_id = apply(rdat_id, 2, as.numeric)
  rdat_id2 = apply(rdat_id2, 2, as.numeric)

  # Row identifiers
  datid = data.frame(as.matrix(1:nrow(dat)))
  dat$id = 1:nrow(dat)
  colnames(datid)[1] = "id"

  # Columns flagging with '+' the proteins to delete (Reverse, contaminants, ...)
  nb_plus = apply(dat, 2, function(x) { sum(x == '+') })
  datdelprot = dat[, which(nb_plus > 0)]
  datdelprot = as.data.frame(datdelprot)
  colnames(datdelprot) = colnames(dat)[which(nb_plus > 0)]

  # 0 = not quantified; log2 transform
  rdat[rdat == 0] = NA
  rdat = log2(rdat)
  rdata = data.frame(rdat)

  list(
    data_init = dat,
    data_fin = rdata,
    data_fin_id = rdat_id,
    data_fin_idd = rdat_id2,
    datdelprot = datdelprot,
    data_id = datid,
    chemin = name,
    sizee = size,
    spectronaut = spectronaut,
    intensity = intensity
  )
}

# Sample columns available for the analysis (to fill cond1_samples / cond2_samples)
list_samples <- function(file, spectronaut = FALSE, intensity = c("Raw", "LFQ")) {
  names(read_input(file, spectronaut, intensity)$data_fin)
}

# Columns flagging proteins to delete with '+'
list_deletion_columns <- function(file, spectronaut = FALSE, intensity = c("Raw", "LFQ")) {
  names(read_input(file, spectronaut, intensity)$datdelprot)
}
