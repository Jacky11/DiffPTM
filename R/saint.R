#################################################
# SAINTexpress (optional): interaction scores for AP-MS experiments
# Condition 1 = control samples, condition 2 = bait samples; the bait name is the
# name of condition 2.
#  - MaxQuant files: preys are identified by their majority protein IDs and protein
#    lengths are read in the "Sequence length" column.
#  - Spectronaut files: preys are identified by the first accession of their protein
#    group, and protein lengths are estimated from PG.MolecularWeight (mean residue
#    mass of 110 Da).
# SAINTexpress is run in the saint/ sub-folder of the output folder, where its
# input and output files are kept.
#################################################

saint_executable <- function(path = NULL) {
  if (!is.null(path) && nzchar(path)) return(normalizePath(path, mustWork = TRUE))
  exe <- if (.Platform$OS.type == "windows") "SAINTexpress-int.exe" else "SAINTexpress-int"
  candidates <- c(file.path(getOption("pairwise.app_dir", "."), "bin", exe), Sys.which(exe))
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(candidates) == 0)
    stop("SAINTexpress not found. Put ", exe, " in the bin/ folder (see bin/README.md) or set 'saint_exe'.")
  normalizePath(candidates[1])
}

run_saint <- function(ds, ldat2, vars1, vars2, p, out_dir, log = function(msg) NULL) {
  exe <- saint_executable(p$saint_exe)
  saint_dir <- file.path(out_dir, "saint")
  dir.create(saint_dir, showWarnings = FALSE)
  c1 <- p$cond1_name; c2 <- p$cond2_name

  dataSaint = ds$data_init
  datSa = ldat2
  if (p$intensity == "LFQ") {
    datSa[, grep("LFQ.intensity.", colnames(datSa))] = 2^datSa[, grep("LFQ.intensity.", colnames(datSa))]
  } else {
    datSa[, grep("Intensity.", colnames(datSa))] = 2^datSa[, grep("Intensity.", colnames(datSa))]
  }

  spectro = isTRUE(p$spectronaut)
  len_col = grep("equence.len", x = colnames(dataSaint))
  if (length(len_col) == 0) {
    mw_col = grep("MolecularWeight", x = colnames(dataSaint))
    if (length(mw_col) == 0)
      stop("SAINT needs the protein lengths: no 'Sequence length' (MaxQuant) or 'PG.MolecularWeight' (Spectronaut) column in the input file.")
    log("SAINT: protein lengths estimated from the molecular weight (PG.MolecularWeight / 110 Da).")
  }
  id_col = if (spectro) "Protein.IDs" else grep("Majority.protein.", x = colnames(dataSaint))

  seq_len = NULL
  maj_prot_Id = NULL
  for (i in 1:nrow(datSa)) {
    if (length(len_col) > 0) {
      req = dataSaint[datSa[i, 1], len_col]
      if (length(req) == 1) {
        seq_len = c(seq_len, as.numeric(req))
      } else {
        seq_len = c(seq_len, as.numeric(unlist(strsplit(as.character(req[1]), ";"))[1]))
      }
    } else {
      mw = convert_number(unlist(strsplit(as.character(dataSaint[datSa[i, 1], mw_col[1]]), ";")))
      mw = mw[!is.na(mw)]
      seq_len = c(seq_len, if (length(mw) > 0) round(mw[1] / 110) else NA)
    }
    maj_prot_Id = c(maj_prot_Id, dataSaint[datSa[i, 1], id_col])
  }
  # proteins without molecular weight: median length
  if (length(len_col) == 0) seq_len[is.na(seq_len)] = round(median(seq_len, na.rm = TRUE))
  # gene names (PreyGene column of SAINTexpress): first gene of the protein group
  first_token = function(v) sapply(strsplit(gsub("\"", "", as.character(v)), ";"), function(x) if (length(x) == 0) NA else trimws(x[1]))
  rows_sa = datSa[, 1]
  if (spectro) {
    gene = first_token(dataSaint$Majority.protein.IDs[rows_sa])
  } else if ("Gene.names" %in% colnames(dataSaint)) {
    gene = first_token(dataSaint$Gene.names[rows_sa])
  } else if ("Fasta.headers" %in% colnames(dataSaint)) {
    gene = sapply(as.character(dataSaint$Fasta.headers[rows_sa]), function(h) {
      m = regmatches(h, regexpr("GN=[^ ;]+", h)); if (length(m) == 0) NA else sub("GN=", "", m)
    }, USE.NAMES = FALSE)
  } else {
    gene = rep(NA, length(rows_sa))
  }
  dataSa = data.frame(datSa, seq = seq_len, prot_id = maj_prot_Id)

  if (p$intensity == "LFQ") {
    dataIntensity = dataSa[, grep("LFQ.intensity.", colnames(dataSa))]
    short.name = unlist(strsplit(colnames(dataIntensity), split = "LFQ.intensity."))
  } else {
    dataIntensity = dataSa[, grep("Intensity.", colnames(dataSa))]
    short.name = unlist(strsplit(colnames(dataIntensity), split = "Intensity."))
  }
  prey = cbind(dataSa$prot_id, dataSa$seq)
  short.name = short.name[short.name != ""]

  # bait file
  databait = data.frame(nom = NULL, type = NULL, test_control = NULL)
  bait_name = gsub("[^A-Za-z0-9_.-]", "_", c2)
  for (i in 1:length(short.name)) {
    is_bait = colnames(dataIntensity)[i] %in% vars2
    databait[i, "test_control"] = if (is_bait) "T" else "C"
    databait[i, "type"] = if (is_bait) bait_name else "Ctrl"
    databait[i, "nom"] = short.name[i]
  }
  databait = databait[, c(3, 2, 1)]
  bait_file  <- file.path(saint_dir, paste0("bait_", c1, "_vs_", c2, ".dat"))
  prey_file  <- file.path(saint_dir, paste0("prey_", c1, "_vs_", c2, ".dat"))
  inter_file <- file.path(saint_dir, paste0("inter_", c1, "_vs_", c2, ".dat"))
  list_file  <- file.path(saint_dir, paste0("list_", c1, "_vs_", c2, ".txt"))
  write.table(databait, bait_file, sep = "\t", row.names = F, col.names = F, quote = F)

  # prey file
  colnames(dataIntensity) = short.name
  baits = databait[, "type"]
  intensities = as.vector(as.matrix(dataIntensity))
  replicate = NULL
  bait = NULL
  for (i in 1:ncol(dataIntensity)) {
    replicate = c(replicate, rep(short.name[i], nrow(dataIntensity)))
    bait = c(bait, rep(baits[i], nrow(dataIntensity)))
  }
  name_mod = rep(NA, nrow(prey))
  for (i in 1:nrow(prey)) {
    aa = gsub(pattern = "\"", replacement = "", x = unlist(strsplit(prey[i, 1], split = ";")))
    name_mod[i] = aa[1]
  }
  # the gene names are not given to SAINTexpress (preys sharing a gene name would be
  # merged); they are added to the PreyGene column of the results
  gene[is.na(gene) | gene == ""] = name_mod[is.na(gene) | gene == ""]
  prey2 = cbind(name_mod, dataSa$seq)
  prey22 = na.omit(prey2)
  prey222 = as.numeric(prey22[, 2])
  prey22 = data.frame(prey22[, 1], prey222)
  write.table(prey22, prey_file, sep = "\t", row.names = F, col.names = F, quote = F)

  # interaction file
  idd = rep(name_mod, ncol(dataIntensity))
  interactions = data.frame(replicate, bait, idd, intensities)
  interactions = interactions[(interactions[, (ncol(interactions))] != 0), ]
  data_fin = na.omit(interactions)
  write.table(data_fin, inter_file, sep = "\t", row.names = F, col.names = F, quote = F)

  # SAINTexpress writes list.txt in the working directory
  old <- setwd(saint_dir); on.exit(setwd(old), add = TRUE)
  status <- system2(exe, basename(c(inter_file, prey_file, bait_file)), stdout = "SAINTexpress.log", stderr = "SAINTexpress.log")
  setwd(old)
  if (!file.exists(file.path(saint_dir, "list.txt")))
    stop("SAINTexpress failed (exit status ", status, "), see ", file.path(saint_dir, "SAINTexpress.log"))
  file.rename(file.path(saint_dir, "list.txt"), list_file)

  dataSAINTFINAL = read.table(list_file, sep = "\t", header = TRUE, fill = TRUE, comment.char = "")
  if ("PreyGene" %in% colnames(dataSAINTFINAL)) {
    g = gene[match(dataSAINTFINAL$Prey, name_mod)]
    dataSAINTFINAL$PreyGene = ifelse(is.na(g), dataSAINTFINAL$PreyGene, g)
  }

  # fold-change between baits and controls (on the filtered, normalised data)
  dataAvgP = dataSAINTFINAL
  dataFold = ldat2
  for (i in 1:ncol(dataFold)) { dataFold[, i] = as.numeric(dataFold[, i]) }
  dataFold[dataFold == 0] = NA

  # bait samples = condition 2, control samples = condition 1
  test = dataFold[, vars2, drop = FALSE]
  control = dataFold[, vars1, drop = FALSE]
  control = apply(control, 1, mean, na.rm = T)
  test = apply(test, 1, mean, na.rm = T)
  foldchange = test - control

  dataAvgP = data.frame(dataAvgP["Prey"], dataAvgP["AvgP"])
  if (spectro) {
    dataprot = dataSaint[, c("Protein.IDs", "Fasta.headers")]
  } else {
    dataprot = dataSaint[, grepl("Majority.protein.|Fasta.headers", colnames(dataSaint))]
  }
  dataprot = dataprot[dataFold[, 1], ]
  dataprot = data.frame(dataprot, foldchange)
  first_id = sapply(strsplit(gsub("\"", "", as.character(dataprot[, 1])), ";"), function(v) v[1])

  dataFinal = data.frame(nameProteins = NULL, FoldChange = NULL, AvgP = NULL, Gene = NULL)
  for (i in 1:nrow(dataAvgP)) {
    x = which(first_id == dataAvgP[i, 1])
    if (length(x) != 0) {
      dataFinal = rbind(dataFinal, data.frame(nameProteins = dataAvgP[i, 1], FoldChange = dataprot[x, 3], AvgP = as.numeric(dataAvgP[i, 2]), Gene = dataprot[x, 2]))
    }
  }

  # NB: the selection uses FoldChange > 1 (log2), while the plot uses log2(fold_change)
  interestingData = subset(dataFinal, FoldChange > 1 & AvgP > p$saint_avgp)
  interestingData[, 4] = gsub("\"", "", interestingData[, 4])
  dataFoldChangemax = head(dplyr::arrange(interestingData, dplyr::desc(FoldChange)), 10)
  dataAvgpmax = head(dplyr::arrange(interestingData, dplyr::desc(AvgP)), 10)
  dataFoldChangemax$Gene = genesearch(dataFoldChangemax)
  dataAvgpmax$Gene = genesearch(dataAvgpmax)
  labelsgraph = rbind(dataAvgpmax, dataFoldChangemax)
  labelsgraph = labelsgraph[!duplicated(labelsgraph, fromLast = T), ]
  if (length(labelsgraph$Gene) != 0) {
    for (i in 1:length(labelsgraph$Gene)) {
      if (is.null(labelsgraph[i, 4][[1]])) {
        labelsgraph[i, 4] = substr(labelsgraph[i, 1], 1, ifelse(nchar(labelsgraph[i, 1]) > 13, 13, nchar(labelsgraph[i, 1])))
      }
    }
  }

  plot = ggplot(dataFinal, aes(y = AvgP, x = FoldChange)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0, 1.2)) +
    scale_x_continuous(breaks = seq(from = floor(min(dataFinal["FoldChange"], na.rm = T)), ceiling(max(dataFinal["FoldChange"], na.rm = T)), by = 1),
                       limits = c(min(dataFinal["FoldChange"], na.rm = T), max(dataFinal["FoldChange"] + 1, na.rm = T))) +
    geom_hline(yintercept = p$saint_avgp, linetype = "dashed", color = "red") +
    ggtitle("Only proteins quantified in both conditions") +
    geom_vline(xintercept = log2(p$fold_change), linetype = "dashed", color = "red") +
    geom_point(na.rm = T, color = ifelse(dataFinal$FoldChange > log2(p$fold_change) & dataFinal$AvgP > p$saint_avgp, "black", "grey")) +
    ylab("AvgP") + xlab(paste0("log2([Mean in ", c2, "] / [Mean in ", c1, "])")) +
    geom_label_repel(data = labelsgraph, mapping = aes(label = Gene), color = "darkgreen", na.rm = T, max.overlaps = 100,
                     box.padding = 1, alpha = 0.7, max.time = 2, max.iter = 20000, size = 2.5) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(color = "grey"), panel.grid.minor = element_line(color = "lightgrey"))

  list(dataSAINTFINAL = dataSAINTFINAL, databait = databait, dataFinal = dataFinal,
       interestingData = interestingData, plot = plot,
       files = c(bait_file, prey_file, inter_file, list_file, file.path(saint_dir, "SAINTexpress.log")))
}
