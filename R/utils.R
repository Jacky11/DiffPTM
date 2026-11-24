#################################################
# Small helper functions
#################################################

# Find the common prefix of a vector of strings
find_common_prefix <- function(strings) {
  if (length(strings) == 0) return("")
  prefix <- strings[[1]]
  for (string in strings[-1]) {
    while (!startsWith(string, prefix)) {
      prefix <- substr(prefix, 1, nchar(prefix) - 1)
      if (nchar(prefix) == 0) return("")
    }
  }
  return(prefix)
}

# Reverse a string
reverse_string <- function(s) {
  return(sapply(lapply(strsplit(s, NULL), rev), paste, collapse = ""))
}

# Convert strings to numbers (accepts both "," and "." as decimal mark)
convert_number <- function(x) {
  x <- as.character(x)
  x <- gsub(pattern = ",", replacement = ".", x = x, fixed = TRUE)
  x <- as.numeric(x)
  return(x)
}

# Maximum value in a semicolon-separated list
find_max <- function(col) {
  values <- unlist(strsplit(as.character(col), ";"))
  values <- convert_number(values)
  max_value <- max(values, na.rm = TRUE)
  return(max_value)
}

contains_semicolon <- function(cell) {
  grepl(";", cell)
}

# Extract the gene name (GN=...) from a FASTA header
find_gene_fasta <- function(x1) {
  x1 = as.character(x1)
  string_split = strsplit(x1, " ")
  result = sapply(string_split, function(x) {
    index = grep("GN=", x)
    if (length(index) > 0) {
      return(substr(x[index], 4, nchar(x[index])))
    }
  })
  result = unique(result)
  result = result[!is.na(result)]
  if (length(result) > 1) {
    for (i in 2:length(result)) {
      result[1] = paste0(result[1], ",", result[i])
      result[i] = NA
    }
  }
  result = result[!is.na(result)]
  return(result)
}

# Gene names from the "Gene" column of the SAINT tables
genesearch <- function(data) {
  return(sapply(data$Gene, function(x1) {
    string_split = strsplit(x1, " ")
    result = sapply(string_split, function(x) {
      index = grep("GN=", x)
      if (length(index) > 0) {
        return(substr(x[index], 4, nchar(x[index])))
      }
    })
    result = unique(result)
    result = result[!is.na(result)]
    if (length(result) > 1) {
      for (i in 2:length(result)) {
        result[1] = paste0(result[1], ",", result[i])
        result[i] = NA
      }
    }
    result = result[!is.na(result)]
    return(result)
  }))
}

# Uniprot accession / protein name / gene name of a set of protein groups.
#   protein_ids   : "Protein.IDs" column (e.g. "sp|P12345|NAME_HUMAN;...")
#   fasta_headers : "Fasta.headers" column
annotate_proteins <- function(protein_ids, fasta_headers) {
  access = as.character(fasta_headers)
  protname = access
  genename = access
  access = protein_ids
  if (length(access) > 0) {
    for (i in 1:length(access)) {
      if (!is.na(access[i])) {
        car = unlist(strsplit(access[i], "[|]"))
        if (length(car) > 1) {
          if (length(car) > 2) {
            colnu = unlist(strsplit(car[3], "[ ]"))
            protname[i] = colnu[1]
            colnn = colnu[grep("GN=", colnu)]
            if (length(colnn) != 0) {
              genename[i] = unlist(strsplit(colnn, "[=]"))[2]
            }
          }
        } else {
          genename[i] = find_gene_fasta(genename[i])
          if (length(genename[i]) == 0) {
            if (length(unlist(strsplit(unlist(strsplit(car[1], "[ ]")), "[>]"))) > 1) {
              access[i] = unlist(strsplit(unlist(strsplit(car[1], "[ ]")), "[>]"))[2]
            }
          }
        }
      }
    }
  }
  list(access = access, protname = protname, genename = genename)
}

# Add the 3 annotation columns in front of a table (as done for every sheet)
add_annotation_columns <- function(tab, ann) {
  tab = data.frame(tab)
  tab$access = ann$access
  tab$protname = ann$protname
  tab$genename = ann$genename
  tab = tab[, c((ncol(tab) - 2):ncol(tab), 1:(ncol(tab) - 3))]
  colnames(tab)[1:3] = c("Uniprot accession", "Protein name", "Gene name")
  tab
}

# Add the columns expected by the output tables when they are missing
# (e.g. Spectronaut files have no "Peptides" or "Mol..weight..kDa." column)
fill_missing_columns <- function(d) {
  if (is.null(d$Fasta.headers)) {
    d = cbind(d, rep(NA, nrow(d)))
    colnames(d)[length(d)] = "Fasta.headers"
  }
  if (is.null(d$Majority.protein.IDs)) {
    d = cbind(d, rep(NA, nrow(d)))
    colnames(d)[length(d)] = "Majority.protein.IDs"
    if (!is.null(d$Protein)) {
      d$Majority.protein.IDs = d$Protein
    }
  }
  if (is.null(d$Peptides)) {
    d = cbind(d, rep(NA, nrow(d)))
    colnames(d)[length(d)] = "Peptides"
  }
  if (is.null(d$Mol..weight..kDa.)) {
    d = cbind(d, rep(NA, nrow(d)))
    colnames(d)[length(d)] = "Mol..weight..kDa."
  }
  d
}

# Short condition name used in Excel sheet names (31 characters max)
short_name <- function(name) {
  if (nchar(name) > 6) paste0(substr(name, 1, 3), substr(name, nchar(name) - 3, nchar(name))) else name
}

# Evaluate an expression without printing what it prints (some helper functions print
# intermediate results)
quiet <- function(expr) {
  utils::capture.output(res <- suppressMessages(expr))
  res
}

# Title of a slide
slide_title <- function(doc, text, size = 28) {
  officer::ph_with(x = doc, location = officer::ph_location_type(type = "title"),
                   value = officer::block_list(officer::fpar(officer::ftext(text, officer::fp_text(font.size = size)))))
}
