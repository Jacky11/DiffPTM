#################################################
# Differential analysis between two conditions
#
# run_pairwise() runs the whole analysis (filters, normalisation, imputation,
# tests, FDR control, figures and Excel sheets). It is called by the Shiny
# application (server.R) and by scripts/run_analysis.R.
#
#   p        : list of parameters (see default_params() in R/params.R)
#   out_dir  : folder where the result files are written
#   data     : optional, result of read_input() (avoids reading the file twice in Shiny)
#   progress : optional function(message) called at each step (Shiny progress bar)
#   log      : function(message) used for log messages
#
# Returns (invisibly) a list with the paths of the files written and the main tables.
#################################################

run_pairwise <- function(p, out_dir = ".", data = NULL,
                         progress = function(msg) NULL,
                         log = function(msg) message(msg)) {

  p <- check_params(p)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- function(f) file.path(out_dir, f)
  step <- function(msg) { log(msg); progress(msg) }
  if (!is.null(p$seed)) set.seed(p$seed)

  Ts = proc.time()

  #################################
  # Data
  #################################
  if (is.null(data) && (is.null(p$input_file) || !file.exists(p$input_file)))
    stop("Input file not found: ", p$input_file)
  ds <- if (is.null(data)) read_input(p$input_file, p$spectronaut, p$intensity) else data
  if (!identical(ds$intensity, p$intensity) || !identical(ds$spectronaut, isTRUE(p$spectronaut)))
    stop("'data' was read with other 'intensity' / 'spectronaut' settings than 'p'.")
  available <- names(ds$data_fin)
  miss <- setdiff(c(p$cond1_samples, p$cond2_samples), available)
  if (length(miss) > 0)
    stop("Unknown sample column(s): ", paste(miss, collapse = ", "),
         "\nAvailable columns: ", paste(available, collapse = ", "))

  del_cols <- if (is.null(p$delete_columns)) names(ds$datdelprot) else p$delete_columns
  miss <- setdiff(del_cols, names(ds$datdelprot))
  if (length(miss) > 0) stop("Unknown deletion column(s): ", paste(miss, collapse = ", "))

  vars1 <- p$cond1_samples
  vars2 <- p$cond2_samples
  nbc1 = length(vars1)
  nbc2 = length(vars2)
  FCT = log2(p$fold_change)
  FDRT = p$fdr
  th <- filter_thresholds(p, nbc1, nbc2)
  ftxt <- filter_methods_text(p, th, p$cond1_name, p$cond2_name)
  ibaq <- isTRUE(p$ibaq)
  editable_file <- out(paste0("Editable_plots_", p$cond2_name, "_vs_", p$cond1_name, ".pptx"))
  editable <- function(x) {
    if (isTRUE(p$editable_pptx)) export::graph2ppt(x, file = editable_file, append = TRUE)
    invisible(NULL)
  }

  #################################
  # PowerPoint: title and description of the pipeline
  #################################
  step("Writing data analysis pipeline")
  doc <- officer::read_pptx()
  doc <- officer::add_slide(doc, layout = "Section Header", master = "Office Theme")
  doc = officer::ph_with(x = doc, location = officer::ph_location_type(type = "title"), value = "Differential Analysis")

  doc <- officer::add_slide(doc, layout = "Title and Content", master = "Office Theme")
  doc = officer::ph_with(x = doc, location = officer::ph_location_type(type = "title"), value = "Data analysis pipeline")

  point0 = paste("Data used:", ds$chemin, "(", ds$sizee, " bytes)")
  point1 = paste0("1. Quantified values: ", ftxt$definition, " ", if (p$intensity == "LFQ") "LFQ" else "Raw",
                  " intensities are log2-transformed; zero intensities are considered as missing values.")
  del_for_text <- if (is.null(p$delete_columns)) names(ds$datdelprot) else p$delete_columns
  point2 = if (length(del_for_text) > 0) {
    paste0("2. Proteins flagged with '+' in the column(s) ", paste(gsub(".", " ", del_for_text, fixed = TRUE), collapse = ", "),
           " of the input file (e.g. reverse sequences, potential contaminants, proteins only identified by site) are removed.")
  } else {
    "2. No column of the input file flags proteins to remove (reverse sequences, contaminants)."
  }
  point3 = paste("3. Filter 1 (number of quantified values):", ftxt$tested)
  if (p$normalization == "Median") {
    point4 = "4. Normalisation: within each condition, the log2 intensities of each sample are centred on the mean of the medians of the samples of the condition. The medians are computed on the proteins kept by the filters (tested and present/absent proteins) that have no missing value in the condition."
  } else {
    point4 = "4. In each condition the intensity values are not normalized."
  }
  point5 = paste("5. Filter 2 (present/absent proteins):", ftxt$present, ftxt$removed)
  if (p$imputation == "MLE") {
    point6 = "6. Imputation: the missing values of the tested proteins are imputed condition by condition with the impute.mle function of the R package imp4p."
  } else if (p$imputation == "SLSA") {
    point6 = "6. Imputation: the missing values of the tested proteins are imputed with the impute.slsa function of the R package imp4p."
  } else {
    point6 = "6. The missing values are not imputed."
  }
  if (p$test == "LIMMA") {
    point7 = "7. Differential analysis is performed. First, a LIMMA t-test is used and the corresponding p-values are next adjusted using an adaptive Benjamini-Hochberg correction."
  } else {
    point7 = "7. Differential analysis is performed. First, a paired LIMMA t-test is used and the corresponding p-values are next adjusted using an adaptive Benjamini-Hochberg correction."
  }
  point8 = paste("The FDR threshold is fixed to", FDRT, "and the proportion of differentially abundant proteins is estimated with the", p$pi0_method,
                 "method thanks to the R package cp4p. The proteins are considered significantly regulated if their adjusted p-value is below ", FDRT,
                 " and if their absolute log2(fold-change) is superior to", FCT, ".")
  if (isTRUE(p$conditional_pvalues))
    point8 = paste(point8, "The p-values are adjusted among the proteins whose absolute log2(fold-change) is superior to", FCT,
                   "after being made conditional on this selection: each p-value is divided by the p-value the protein would have with an absolute log2(fold-change) equal to", FCT, ".")

  fp_normal <- officer::fp_text(font.size = 14)
  pars <- officer::block_list(
    officer::fpar(officer::ftext(point0, fp_normal)), officer::fpar(officer::ftext(point1, fp_normal)),
    officer::fpar(officer::ftext(point2, fp_normal)), officer::fpar(officer::ftext(point3, fp_normal)),
    officer::fpar(officer::ftext(point4, fp_normal)), officer::fpar(officer::ftext(point5, fp_normal)),
    officer::fpar(officer::ftext(point6, fp_normal)), officer::fpar(officer::ftext(point7, fp_normal)),
    officer::fpar(officer::ftext(point8, fp_normal)))
  doc = officer::ph_with(x = doc, location = officer::ph_location(left = 0.6, top = 2, height = 5, width = 9), value = pars)

  #################################
  # Proteins to delete (flagged by '+') and parameters
  #################################
  dat = ds$data_init
  idc = ds$data_id[, 1]

  if (length(del_cols) > 0) {
    sdelprot = apply(ds$datdelprot[, del_cols, drop = FALSE], 1, function(x) sum(x == '+'))
  } else {
    sdelprot = 0
  }
  delprot = (sdelprot > 0)


  #################################
  # Sample names
  #################################
  strip_prefix <- function(v) {
    v = strsplit(x = v, split = "[.]")
    res = NULL
    first = if (p$intensity == "Raw") 2 else 3
    for (i in 1:length(v)) { res = c(res, paste(v[[i]][first:length(v[[i]])], collapse = ".")) }
    res
  }
  name_c1 = strip_prefix(vars1)
  name_c2 = strip_prefix(vars2)

  pepnamec1 = paste("Peptides.", name_c1, sep = "")
  pepnamec2 = paste("Peptides.", name_c2, sep = "")
  upepnamec1 = paste("Unique.peptides.", name_c1, sep = "")
  upepnamec2 = paste("Unique.peptides.", name_c2, sep = "")

  # columns needed for each selected sample
  miss <- c(setdiff(c(pepnamec1, pepnamec2), colnames(ds$data_fin_id)),
            setdiff(c(upepnamec1, upepnamec2), colnames(ds$data_fin_idd)))
  if (length(miss) > 0)
    stop("Column(s) not found in the input file: ", paste(miss, collapse = ", "),
         ". Each selected sample needs 'Peptides <sample>' and 'Unique peptides <sample>' columns.")
  if (ibaq) {
    miss <- setdiff(paste0("iBAQ.", c(name_c1, name_c2)), colnames(dat))
    if (length(miss) > 0) {
      log(paste0("Warning: iBAQ column(s) not found (", paste(miss, collapse = ", "),
                 "): the analysis is run without iBAQ values."))
      ibaq <- FALSE
      p$ibaq <- FALSE
    }
  }
  parameters = parameters_table(p, del_cols)

  #################################
  # Data matrices; a value is quantified if it has enough unique peptides
  #################################
  rdata = as.matrix(cbind(idc, as.matrix(ds$data_fin[, vars1]), as.matrix(ds$data_fin[, vars2])))

  if (sum(delprot) == 0) { rdat_id1 = as.matrix(ds$data_fin_id[, pepnamec1]) } else { rdat_id1 = as.matrix(ds$data_fin_id[!delprot, pepnamec1]) }
  if (sum(delprot) == 0) { rdat_id2 = as.matrix(ds$data_fin_id[, pepnamec2]) } else { rdat_id2 = as.matrix(ds$data_fin_id[!delprot, pepnamec2]) }

  rdat_Unique = as.matrix(cbind(as.matrix(ds$data_fin_idd[, upepnamec1]), as.matrix(ds$data_fin_idd[, upepnamec2])))
  rdata[, 2:ncol(rdata)][which(rdat_Unique < p$min_unique_peptides)] = NA

  #################################
  # Venn diagrams / UpSet plots
  #################################
  step("Venn diagrams")
  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  n0 = nrow(rdata)
  rdata_rr = rdata[!delprot, ]
  doc = slide_title(doc, "Number of identified proteins by replicates and conditions after deleting reverse, contaminants and only identified by site proteins")

  sample_part <- function(m) {
    b = NULL
    for (i in 1:ncol(m)) {
      a = unlist(strsplit(as.character(colnames(m)[i]), '.', fixed = TRUE))
      b = c(b, paste(a[2:length(a)], collapse = "."))
    }
    b
  }
  colnames(rdat_id1) = sample_part(rdat_id1)
  v1 <- quiet(makeVenn(rdat_id1, nintersects = 25))
  editable(v1)
  doc <- officer::ph_with(doc, value = v1, location = officer::ph_location(left = 0, top = 2, height = 5, width = 5))

  colnames(rdat_id2) = sample_part(rdat_id2)
  v2 <- quiet(makeVenn(rdat_id2, nintersects = 25))
  editable(v2)
  doc <- officer::ph_with(doc, value = v2, location = officer::ph_location(left = 5, top = 2, height = 5, width = 5))

  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  doc = slide_title(doc, "Number of identified proteins by replicates and conditions after deleting reverse, contaminants and only identified by site proteins")
  rdat_all = cbind(rdat_id1, rdat_id2)
  vall <- quiet(makeVenn(rdat_all, nintersects = 25))
  doc <- officer::ph_with(doc, value = vall, location = officer::ph_location(left = 0, top = 2, height = 6, width = 10))
  editable(vall)

  n3 = nrow(rdata_rr)
  ldat = rdata_rr
  colnames(ldat)[1] = "id"

  #################################
  # Missing values and distribution of intensities
  #################################
  step("Missing values and distributions")
  coln = colnames(ldat[, 2:ncol(ldat)])

  conditions = as.factor(c(rep(p$cond1_name, length(vars1)), rep(p$cond2_name, length(vars2))))
  conditions = factor(conditions, levels(conditions)[c(which(levels(conditions) == p$cond1_name), which(levels(conditions) == p$cond2_name))])

  nb_NA_repl = apply(ldat[, 2:ncol(ldat)], 2, function(x) sum(is.na(x)))
  nb_NA_prot = apply(ldat[, 2:ncol(ldat)], 1, function(x) sum(is.na(x)))

  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  doc = slide_title(doc, "Missing values and distribution of quantified values by sample")

  df = data.frame(nb_NA_repl, coln, conditions)
  colnames(df) = c("NA_per_sample", "Name", "Conditions")
  ggp = ggplot(data = df, aes(x = Name, y = NA_per_sample, fill = Conditions)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = NA_per_sample), hjust = 1.6, color = "black", size = 3.5, angle = 90) +
    labs(x = "", y = paste0("Missing values per sample\n (", nrow(ldat), " proteins identified among samples)")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")
  doc <- officer::ph_with(doc, value = ggp, location = officer::ph_location_left())
  editable(ggp)

  dataplot = ldat[which(nb_NA_prot == 0), 2:ncol(ldat)]
  dfd = stack(data.frame(dataplot))
  dfd = cbind(dfd, conditions)
  colnames(dfd) = c("log2_intensities", "Samples", "Conditions")
  ggp <- ggplot(dfd, aes(x = log2_intensities, y = Samples, fill = after_stat(x))) +
    stat_density_ridges(geom = "density_ridges_gradient", scale = 3, rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 0.5, alpha = 0.7) +
    scale_fill_viridis_c(name = "log2(int.)", option = "C", direction = -1) +
    labs(title = "Distribution of intensities for proteins\n without missing values")
  doc <- officer::ph_with(doc, value = ggp, location = officer::ph_location_right())
  editable(ggp)

  # same slide without colours (version "ALL")
  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  doc = slide_title(doc, "Missing values and distribution of quantified values by sample")
  ggp2 = ggplot(data = df, aes(x = Name, y = NA_per_sample)) +
    geom_bar(stat = "identity", fill = "lightblue") +
    geom_text(aes(label = NA_per_sample), hjust = 1.6, color = "black", size = 3.5, angle = 90) +
    labs(x = "", y = paste0("Missing values per sample\n (", nrow(ldat), " proteins identified among samples)")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
  doc <- officer::ph_with(doc, value = ggp2, location = officer::ph_location_left())
  editable(ggp2)
  doc <- officer::ph_with(doc, value = ggp, location = officer::ph_location_right())
  editable(ggp)

  #################################
  # Clustering of samples on quantified values (Pearson correlation, Ward.D2)
  #################################
  step("Correlation matrices")
  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  doc = slide_title(doc, "Clustering of samples using quantified values")

  z <- cor(ldat[, 2:ncol(ldat)], use = "pairwise.complete.obs")
  d <- corr_heatmap(z, gradientRate = 2.5)
  matrice_corr = z
  doc <- officer::ph_with(doc, value = d, location = officer::ph_location(left = 0, top = 1.2, height = 5, width = 6))
  editable(d)

  hv <- hclust(as.dist(1 - z), method = "ward.D2")
  ggd <- sample_dendrogram(hv, coln, conditions, "Ward's method with a \ncorrelation-based distance")
  doc <- officer::ph_with(doc, value = ggd, location = officer::ph_location(left = 6, top = 2, height = 5, width = 4))
  editable(ggd)
  doc <- note_block(doc, "Pearson correlation coefficient between each pair of samples is computed using all complete pairs of observations in these samples.",
                    paste0(c("Min. corr. value = ", min(z, na.rm = T)), sep = ""))

  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  doc = slide_title(doc, "Clustering of samples using quantified values")
  ggd <- sample_dendrogram(hv, coln, conditions, "Ward's method with a \ncorrelation-based distance", colour = FALSE)
  doc <- officer::ph_with(doc, value = ggd, location = officer::ph_location(left = 6, top = 2, height = 5, width = 4))
  editable(ggd)

  #################################
  # Clustering of samples on missing values (Jaccard index)
  #################################
  step("Clustering of samples using missing values")
  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  doc = slide_title(doc, "Clustering of samples using missing values")

  df2 = data.frame(is.na(ldat[, 2:ncol(ldat)]))
  colnames(df2) = coln
  d = jaccard_missing(df2, 1)
  d[is.nan(d)] = 1
  z = d
  dq <- corr_heatmap(z, gradientRate = 1.2)
  doc <- officer::ph_with(doc, value = dq, location = officer::ph_location(left = 0, top = 1.2, height = 5, width = 6))
  editable(dq)

  hv <- hclust(as.dist(1 - d), method = "ward.D2")
  ggd <- sample_dendrogram(hv, coln, conditions, "Ward's method with a \nJaccard index-based distance")
  doc <- officer::ph_with(doc, value = ggd, location = officer::ph_location(left = 6, top = 2, height = 5, width = 4))
  editable(ggd)
  doc <- note_block(doc, "Jaccard index between each pair of samples is computed by replacing each missing value by 1 and each observed value by 0.",
                    paste0(c("Min. Jaccard index = ", min(z, na.rm = T)), sep = ""))

  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  doc = slide_title(doc, "Clustering of samples using missing values")
  ggd <- sample_dendrogram(hv, coln, conditions, "Ward's method with a \nJaccard index-based distance", colour = FALSE)
  doc <- officer::ph_with(doc, value = ggd, location = officer::ph_location(left = 6, top = 2, height = 5, width = 4))
  editable(ggd)

  #################################
  # Sample overlap and clustering robustness (optional): Jaccard / overlap indices
  # on protein presence and pvclust (AU/BP) on the same tree as the correlation
  # dendrogram above
  #################################
  overlap_files <- character(0)
  overlap_results <- NULL
  if (isTRUE(p$overlap)) {
    step("Jaccard indices and pvclust (AU/BP)")
    x_rev   <- ldat[, 2:ncol(ldat), drop = FALSE]
    lab_rev <- clean_sample_names(colnames(x_rev))
    grp_rev <- make_groups(lab_rev, p$overlap_groups)
    log(paste("Groups:", paste(levels(grp_rev), table(grp_rev), sep = "=", collapse = ", ")))
    ncores  <- as.integer(p$n_cores)
    rev_res <- run_overlap_analyses(
      x_rev, grp_rev,
      prefix   = ifelse(nzchar(p$overlap_prefix), p$overlap_prefix, "Overlap"),
      doc      = doc,
      nboot    = p$overlap_nboot,
      nperm    = p$overlap_nperm,
      parallel = if (ncores > 1) ncores else FALSE,
      alpha    = 0.95,
      seed     = 123,
      outdir   = out_dir,
      modif    = isTRUE(p$editable_pptx),
      modif_file = editable_file)
    doc <- rev_res$doc
    overlap_files <- rev_res$files
    overlap_results <- rev_res$results
  }

  #################################
  # PCA on proteins without missing values
  #################################
  step("PCA between samples without missing values")
  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  pca_title <- paste0("PCA between samples from ", nrow(na.omit(ldat[, 2:ncol(ldat)])), " proteins without m.v.")
  doc = slide_title(doc, pca_title)

  acp_data = prcomp(na.omit(ldat[, 2:ncol(ldat)]), center = TRUE, scale. = TRUE)
  theta <- seq(0, 2 * pi, length.out = 100)
  circle <- data.frame(x = cos(theta), y = sin(theta))

  b = rep(0, nrow(acp_data$rotation))
  for (i in 1:nrow(acp_data$rotation)) {
    a = unlist(strsplit(row.names(acp_data$rotation)[i], '.', fixed = TRUE))
    if (p$intensity == "LFQ") { b[i] = paste0(a[3:length(a)], collapse = "_") }
    if (p$intensity == "Raw") { b[i] = paste0(a[2:length(a)], collapse = "_") }
  }
  loadings = data.frame(acp_data$rotation, .names = b)
  aaa = summary(acp_data)

  df = data.frame(NA_per_sample = aaa$importance[2, ],
                  Name = paste("Component", 1:length(aaa$importance[2, ])),
                  Condit = 1:length(aaa$importance[2, ]))
  df$Name = factor(df$Name, levels = paste("Component", 1:length(aaa$importance[2, ])))
  colnames(df) = c("Proportion_of_explained_variance", "Component", "Conditions")
  ggp = ggplot(data = df, aes(x = Component, y = Proportion_of_explained_variance)) +
    geom_bar(stat = "identity", color = "blue", fill = "white") +
    geom_text(aes(label = Proportion_of_explained_variance), hjust = -0.2, color = "black", size = 3.5, angle = 90) +
    ylim(c(0, 1)) +
    labs(x = "", y = paste0("Proportion of explained variance")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")
  doc <- officer::ph_with(doc, value = ggp, location = officer::ph_location_left())
  editable(ggp)

  pca_plot <- function(pc1, pc2, colour) {
    i1 <- as.integer(sub("PC", "", pc1)); i2 <- as.integer(sub("PC", "", pc2))
    pp <- ggplot(circle, aes(x, y)) + geom_path(colour = "white")
    if (colour) {
      pp = pp + geom_text(data = loadings, mapping = aes(x = .data[[pc1]], y = .data[[pc2]], label = .names, colour = as.factor(conditions), cex = 4))
    } else {
      pp = pp + geom_text(data = loadings, mapping = aes(x = .data[[pc1]], y = .data[[pc2]], label = .names, colour = "black", cex = 4))
    }
    pp + coord_fixed(ratio = 1) +
      labs(x = paste(pc1, "(", round(aaa$importance[2, i1] * 100, 2), "%)"),
           y = paste(pc2, "(", round(aaa$importance[2, i2] * 100, 2), "%)")) +
      theme(legend.position = "none")
  }
  pp <- pca_plot("PC1", "PC2", TRUE)
  doc <- officer::ph_with(doc, value = pp, location = officer::ph_location_right())
  editable(pp)

  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  doc = slide_title(doc, pca_title)
  pp <- pca_plot("PC1", "PC2", FALSE)
  doc <- officer::ph_with(doc, value = pp, location = officer::ph_location_right())
  editable(pp)

  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  doc = slide_title(doc, pca_title)
  if (ncol(acp_data$rotation) >= 3) {
    pp <- pca_plot("PC2", "PC3", TRUE)
    doc <- officer::ph_with(doc, value = pp, location = officer::ph_location_left())
    editable(pp)
  }
  if (ncol(acp_data$rotation) >= 4) {
    pp <- pca_plot("PC3", "PC4", TRUE)
    doc <- officer::ph_with(doc, value = pp, location = officer::ph_location_right())
    editable(pp)
  }

  #################################
  # Filter 1: proteins without enough quantified values
  #################################
  step("Filtering proteins without enough values")
  tab = ldat[, 2:ncol(ldat)]
  k = 1
  llc = length(levels(conditions))
  nb_rep = rep(0, llc)
  dat_nbNA = matrix(0, nrow(tab), llc)
  for (i in 1:llc) {
    nb_rep[i] = sum((conditions == levels(conditions)[i]))
    dat_nbNA[, i] = apply(tab[, (k:(k + nb_rep[i] - 1))], 1, function(x) sum(is.na(x)))
    k = k + nb_rep[i]
  }
  # number of quantified values per condition and category of each protein (R/filters.R)
  obs1 = nb_rep[1] - dat_nbNA[, 1]
  obs2 = nb_rep[2] - dat_nbNA[, 2]
  flt = classify_proteins(obs1, obs2, th, p$cond1_name, p$cond2_name)
  nobs = data.frame(id = ldat[, 1], o1 = obs1, o2 = obs2, reason = flt$reason)
  count_cols = c(paste0("Nb quantified values in ", p$cond1_name, " (out of ", nb_rep[1], ")"),
                 paste0("Nb quantified values in ", p$cond2_name, " (out of ", nb_rep[2], ")"))
  add_counts <- function(tab, ids) {
    if (is.null(ids) || nrow(tab) == 0 || length(ids) != nrow(tab)) return(tab)
    j = match(ids, nobs$id)
    tab[[count_cols[1]]] = nobs$o1[j]
    tab[[count_cols[2]]] = nobs$o2[j]
    tab
  }

  liste_prot_1opc = ldat[flt$category == "removed", 1]
  ldat2 = ldat[!ldat[, 1] %in% liste_prot_1opc, ]
  n4 = nrow(ldat2)

  datdel = dat[dat$id %in% liste_prot_1opc, ]
  if (is.null(datdel$Fasta.headers)) {
    datdel = cbind(datdel, rep(NA, nrow(datdel)))
    colnames(datdel)[length(datdel)] = "Fasta.headers"
  }
  ann = annotate_proteins(datdel$Protein.IDs, datdel$Fasta.headers)
  if (is.null(datdel$Majority.protein.IDs)) {
    datdel = cbind(datdel, rep(NA, nrow(datdel)))
    colnames(datdel)[length(datdel)] = "Majority.protein.IDs"
    if (!is.null(datdel$Protein)) { datdel$Majority.protein.IDs = datdel$Protein }
  }
  if (is.null(datdel$Peptides)) {
    datdel = cbind(datdel, rep(NA, nrow(datdel)))
    colnames(datdel)[length(datdel)] = "Peptides"
  }
  if (is.null(datdel$Mol..weight..kDa.)) {
    datdel = cbind(datdel, rep(NA, nrow(datdel)))
    colnames(datdel)[length(datdel)] = "Mol..weight..kDa."
  }
  int_pattern <- if (p$intensity == "LFQ") "LFQ" else "Intensity"
  datdel_fin <- data.frame(datdel$Fasta.headers, datdel$Majority.protein.IDs, datdel$Mol..weight..kDa., as.matrix(datdel[, grep(int_pattern, names(datdel), value = TRUE)]))
  datdel_fin$access = ann$access
  datdel_fin$protname = ann$protname
  datdel_fin$genename = ann$genename
  datdel_fin <- datdel_fin[, c((ncol(datdel_fin) - 2):ncol(datdel_fin), 1:(ncol(datdel_fin) - 3))]
  colnames(datdel_fin)[1:6] = c("Uniprot accession", "Protein name", "Gene name", "Fasta.headers", "Majority.protein.IDs", "Mol.weight.kDa.")
  if (nrow(datdel_fin) > 0) {
    j = match(datdel$id, nobs$id)
    datdel_fin = data.frame(datdel_fin[, 1:3], Reason = nobs$reason[j], check.names = FALSE,
                            setNames(data.frame(nobs$o1[j], nobs$o2[j]), count_cols),
                            datdel_fin[, 4:ncol(datdel_fin)])
  }

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "without_enough_values")
  openxlsx::writeDataTable(wb, "without_enough_values", x = datdel_fin)

  #################################
  # Normalisation: within each condition, centring on the mean of the sample
  # medians (medians computed on proteins without missing value in the condition)
  #################################
  step("Normalisation")
  if (p$normalization == "Median") {
    moy = rep(0, llc)
    k = 1
    j = 1
    for (i in 1:llc) {
      nb_NA_prot = apply(ldat2[, 2:ncol(ldat2)][, (k:(k + nb_rep[i] - 1))], 1, function(x) sum(is.na(x)))
      medianes = apply(ldat2[, 2:ncol(ldat2)][which(nb_NA_prot == 0), (k:(k + nb_rep[i] - 1))], 2, median, na.rm = T)
      moy[i] = mean(medianes)
      while (j < (k + nb_rep[i])) {
        ldat2[, j + 1] = ldat2[, j + 1] - medianes[j - k + 1] + moy[i]
        j = j + 1
      }
      k = k + nb_rep[i]
    }
  }

  #################################
  # Imputation of missing values
  #################################
  step("Imputing missing data")
  if (p$imputation == "None") {
    ldat2i = ldat2[, 2:ncol(ldat2)]
  }
  if (p$imputation == "SLSA") {
    ldat2i = quiet(imp4p::impute.slsa(ldat2[, 2:ncol(ldat2)], conditions = conditions, progress.bar = FALSE))
  }
  if (p$imputation == "MLE") {
    ldat3 = ldat2[, 2:ncol(ldat2)]
    ldat2i = imp4p::impute.mle(ldat3, conditions = conditions)
  }
  label_impute = dat[dat$id %in% ldat2[, 1], ]
  imput_matrix = data.frame(Majority.protein.Ids = label_impute[, grep("Majority.protein.", x = colnames(label_impute))],
                            Fasta.headers = label_impute$Fasta.headers, ldat2i)
  imput_matrix = add_counts(imput_matrix, label_impute$id)

  #################################
  # SAINT (optional)
  #################################
  saint_res <- NULL
  if (isTRUE(p$saint)) {
    step("SAINT analysis")
    saint_res <- run_saint(ds, ldat2, vars1, vars2, p, out_dir, log = log)
    openxlsx::addWorksheet(wb, "SAINT")
    openxlsx::writeDataTable(wb, "SAINT", x = saint_res$dataSAINTFINAL)
    openxlsx::addWorksheet(wb, "Bait")
    openxlsx::writeDataTable(wb, "Bait", x = saint_res$databait)
  }

  #################################
  # Filter 2: proteins present in one condition and absent in the other
  #################################
  step("Filtering present/absent proteins")
  liste_prot_p1a2 = ldat[flt$category == "present1", 1]
  liste_prot_a1p2 = ldat[flt$category == "present2", 1]

  ldat2e = ldat2[!ldat2[, 1] %in% c(liste_prot_p1a2, liste_prot_a1p2), ]
  n5 = nrow(ldat2e)

  name1 = as.character(levels(conditions)[1])
  name2 = as.character(levels(conditions)[2])
  shortname1 = short_name(name1)
  shortname2 = short_name(name2)

  # proteins present in condition 1 and absent in condition 2
  aviBAQp1a2 = NULL
  if (length(liste_prot_p1a2) > 0) {
    datp1a2 = dat[dat$id %in% liste_prot_p1a2, ]
    ann = annotate_proteins(datp1a2$Protein.IDs, datp1a2$Fasta.headers)
    datp1a2 = add_annotation_columns(datp1a2, ann)
    if (ibaq) {
      ib = mean_log2_ibaq(datp1a2, name_c1)
      aviBAQp1a2 = ib
      datp1a2 = cbind(datp1a2[, 1], aviBAQp1a2, datp1a2[, 2:ncol(datp1a2)])
      colnames(datp1a2)[1] = "Uniprot accession"
      colnames(datp1a2)[2] = paste("Average iBAQ ", p$cond1_name, sep = "")
    }
  } else { datp1a2 = data.frame(0) }
  openxlsx::addWorksheet(wb, paste0("PresIn", shortname1, "_AbsIn", shortname2))
  openxlsx::writeDataTable(wb, paste0("PresIn", shortname1, "_AbsIn", shortname2), x = add_counts(datp1a2, datp1a2$id))

  # proteins absent in condition 1 and present in condition 2
  aviBAQa1p2 = NULL
  if (length(liste_prot_a1p2) > 0) {
    data1p2 = dat[dat$id %in% liste_prot_a1p2, ]
    ann = annotate_proteins(data1p2$Protein.IDs, data1p2$Fasta.headers)
    data1p2 = add_annotation_columns(data1p2, ann)
    if (ibaq) {
      aviBAQa1p2 = mean_log2_ibaq(data1p2, name_c2)
      data1p2 = cbind(data1p2[, 1], aviBAQa1p2, data1p2[, 2:ncol(data1p2)])
      colnames(data1p2)[1] = "Uniprot accession"
      colnames(data1p2)[2] = paste("Average iBAQ ", p$cond2_name, sep = "")
    }
  } else { data1p2 = data.frame(0) }
  openxlsx::addWorksheet(wb, paste0("PresIn", shortname2, "_AbsIn", shortname1))
  openxlsx::writeDataTable(wb, paste0("PresIn", shortname2, "_AbsIn", shortname1), x = add_counts(data1p2, data1p2$id))

  #################################
  # Number of proteins along the pipeline; distributions after norm. and imp.
  #################################
  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  doc = slide_title(doc, "Number of selected proteins along the data analysis pipeline")
  df = data.frame(NA_per_sample = c(n0, n3, n4, n5),
                  Name = c("Initially", "Filter wrongly id. prot.", "Filter prot. without enough values", "Filter present/absent prot."),
                  Condit = c(1:4))
  df$Name = factor(df$Name, levels = df$Name)
  colnames(df) = c("nb_of_selected", "Filter", "Conditions")
  ggp = suppressWarnings(ggplot(data = df, aes(x = Filter, y = nb_of_selected)) +
    geom_bar(stat = "identity", color = "blue", fill = "white") +
    geom_text(aes(label = nb_of_selected), hjust = 1.6, color = "black", size = 3.5, angle = 90) +
    labs(x = "", y = paste0("nb of selected")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom"))
  doc <- officer::ph_with(doc, value = ggp, location = officer::ph_location_left())
  editable(ggp)

  dataplot = ldat2i
  dfd = stack(data.frame(dataplot))
  dfd = cbind(dfd, conditions)
  colnames(dfd) = c("log2_intensities", "Samples", "Conditions")
  ggp <- ggplot(dfd, aes(x = log2_intensities, y = Samples, fill = after_stat(x))) +
    stat_density_ridges(geom = "density_ridges_gradient", scale = 3, rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 0.5, alpha = 0.7) +
    scale_fill_viridis_c(name = "log2(int.)", option = "C", direction = -1) +
    labs(title = "Distribution of intensities for proteins\n without missing values (after norm and imp")
  doc <- officer::ph_with(doc, value = ggp, location = officer::ph_location_right())
  editable(ggp)

  #################################
  # Clustering of samples after normalisation and imputation
  #################################
  step("Clustering of samples after normalisation and imputation")
  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  doc = slide_title(doc, "Clustering of samples using quantified values after normalisation and imputation")
  z <- cor(ldat2i, use = "pairwise.complete.obs")
  d <- corr_heatmap(z, gradientRate = 3.5)
  matrice_corri = z
  doc <- officer::ph_with(doc, value = d, location = officer::ph_location(left = 0, top = 1.2, height = 5, width = 6))
  editable(d)

  hv <- hclust(as.dist(1 - z), method = "ward.D2")
  ggd <- sample_dendrogram(hv, coln, conditions, "Ward's method with a \ncorrelation-based distance")
  doc <- officer::ph_with(doc, value = ggd, location = officer::ph_location(left = 6, top = 2, height = 5, width = 4))
  editable(ggd)
  doc <- note_block(doc, "Pearson correlation coefficient between each pair of samples is computed using all complete pairs of observations in these samples.",
                    paste0(c("Min. corr. value = ", min(z, na.rm = T)), sep = ""))

  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  doc = slide_title(doc, "Clustering of samples using quantified values after normalisation and imputation")
  ggd <- sample_dendrogram(hv, coln, conditions, "Ward's method with a \ncorrelation-based distance", colour = FALSE)
  doc <- officer::ph_with(doc, value = ggd, location = officer::ph_location(left = 6, top = 2, height = 5, width = 4))
  editable(ggd)

  #################################
  # Differential analysis
  #################################
  step("Differential analysis")
  moyy = matrix(0, nrow(ldat2i), llc)
  k = 1
  for (i in 1:llc) {
    moyy[, i] = apply(ldat2i[, (k:(k + nb_rep[i] - 1))], 1, mean, na.rm = T)
    k = k + nb_rep[i]
  }
  FC = moyy[, 2] - moyy[, 1]
  # only the tested proteins are considered (present/absent proteins are reported separately)
  FC[!ldat2[, 1] %in% ldat[flt$category == "tested", 1]] = NA

  if (p$test == "paired.LIMMA") {
    res_ana = quiet(limma.test.paired(ldat2i, Conditions = conditions))
    pval = res_ana$p.value[, ncol(res_ana$p.value)]
  }
  if (p$test == "LIMMA") {
    res_ana = quiet(limma.test(ldat2i, Conditions = conditions, RepBio = as.factor(1:length(conditions)), RepTech = as.factor(1:length(conditions))))
    pval = res_ana$p.value[, ncol(res_ana$p.value)]
  }

  dat_filtered_FC_init = dat[dat$id %in% ldat2[, 1], ]
  dat_filtered_FC = dat_filtered_FC_init[which(abs(FC) > FCT), ]

  # p-values of the proteins passing the fold-change threshold are adjusted with
  # the adaptive Benjamini-Hochberg procedure (pi0 estimated with cp4p).
  # Optionally, the p-values are first made conditional on the selection
  # |log2FC| > threshold: p_cond = p / c, where c is the p-value the protein would
  # have with |log2FC| equal to the threshold (given its moderated standard error).
  # Under H0, p_cond is uniform among the selected proteins.
  sel_fc = which(abs(FC) > FCT)
  pvalu = pval[sel_fc]
  FCs = FC[sel_fc]
  if (isTRUE(p$conditional_pvalues)) {
    cutoff = conditional_cutoff(res_ana, FCT)
    pval_cond = pmin(1, pval / cutoff)
    res_abh = quiet(cp4p::adjust.p(pval_cond[sel_fc], pi0.method = p$pi0_method))
  } else {
    res_abh = quiet(cp4p::adjust.p(pvalu, pi0.method = p$pi0_method))
  }
  # raw p-value threshold corresponding to the FDR level (horizontal line of the volcano plot)
  pval_threshold = suppressWarnings(max(pvalu[which(res_abh$adjp[, 2] < FDRT)]))

  # first volcano plot
  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  da_title <- paste("Differential Analysis with log2(Fold-Change) =", floor(FCT * 100) / 100, "and FDR =", floor(FDRT * 1000) / 10, "%.")
  doc = slide_title(doc, da_title)

  dat_prob = dat[dat$id %in% ldat2[, 1], ]
  code = dat_prob$Fasta.headers
  colnn = as.character(code)
  for (i in 1:length(colnn)) {
    car = unlist(strsplit(colnn[i], "[|]"))
    if (length(car) > 2) { colnn[i] = car[2] }
  }
  colnn[is.na(colnn)] = "unknown"
  colnn2 = colnn
  ann = annotate_proteins(dat_prob$Protein.IDs, code)
  genename = ann$genename
  genename[is.na(genename)] = "unknown"

  log2FoldChange = FC
  Uniprot.name = colnn2
  Significant = rep(0, length(FC))
  proteins = data.frame(Significant, Uniprot.name, log2FoldChange, pval, rep(1, length(FC)))
  colnames(proteins) = c("significant", "name", "log2FC", "pvalue", "padj")
  proteins$padj[which(abs(FC) > FCT)] = res_abh$adjp[, 2]
  proteins$Significant <- ifelse(proteins$padj < FDRT, "Significant", "Not Significant")
  a = ggplot(proteins, aes(x = log2FC, y = -log10(pvalue))) +
    geom_point(aes(color = Significant)) +
    geom_vline(aes(xintercept = FCT), color = "grey") +
    geom_vline(aes(xintercept = -FCT), color = "grey") +
    geom_hline(aes(yintercept = -log10(pval_threshold)), color = "grey") +
    scale_color_manual(values = c("grey80", "#79a7fc")) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom") +
    theme(legend.title = element_blank()) +
    xlab(paste("log2([ Mean in", conditions[1 + nb_rep[1]], "]/[ Mean in", conditions[1], "])")) +
    theme(axis.title.y = element_text(size = 9)) +
    theme(axis.title.x = element_text(size = 9))
  proteins$name = genename
  doc <- officer::ph_with(doc, value = a, location = officer::ph_location(left = .5, top = 1.8, height = 4.5, width = 5))
  editable(a)

  # second volcano plot (+ average iBAQ of present/absent proteins)
  step("Volcano plot")
  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  doc = slide_title(doc, da_title)

  datim = proteins[, 2:5]
  colnames(datim) = c("names", "logFC", "P.Value", "AP")
  datim[, 3] = -log10(datim[, 3])
  sel = 1 * (datim[, 4] < FDRT)
  sel[is.na(sel)] = 0
  sel[which((sel == 1) & (datim$logFC < (-FCT)))] = -1
  sel[which((sel == 1) & (datim$logFC < (FCT)))] = 0
  datim = cbind(datim, sel)
  colnames(datim)[5] = "Differential"
  datim$Differential = as.numeric(datim$Differential)
  size_point = rep(2.5, nrow(datim))
  size_point[sel != 0] = 3
  cols <- c("0" = "grey90", "-1" = "chartreuse3", "1" = "dodgerblue", "2" = "indianred1")

  # gene names
  access = as.character(datim$names)
  genename = access
  if (length(access) > 0) {
    for (i in 1:length(access)) {
      if (!is.na(access[i])) {
        car = unlist(strsplit(access[i], "[|]"))
        if (length(car) > 1) {
          if (length(car) > 2) {
            colnu = unlist(strsplit(car[3], "[ ]"))
            colnn = colnu[grep("GN=", colnu)]
            if (length(colnn) != 0) { genename[i] = unlist(strsplit(colnn, "[=]"))[2] }
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
  datim$names = genename

  volflu = ggplot(datim, aes(x = logFC, y = P.Value, fill = as.factor(Differential))) +
    ggtitle(label = " ") +
    scale_fill_manual(values = cols) +
    geom_vline(xintercept = -FCT, linetype = "dashed", color = "grey50", linewidth = 1) +
    geom_vline(xintercept = FCT, linetype = "dashed", color = "grey50", linewidth = 1) +
    geom_hline(yintercept = min(datim$P.Value[sel != 0], na.rm = T), linetype = "dashed", color = "grey50", linewidth = 1) +
    geom_text_repel(data = subset(datim, AP < FDRT), aes(label = names), size = 2, segment.color = 'grey50') +
    geom_point(size = size_point, alpha = 1, na.rm = T, shape = 21, colour = "grey70") +
    theme_bw(base_size = 14) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(face = "bold", colour = "black", size = 12)) +
    theme(axis.text.y = element_text(face = "bold", colour = "black", size = 12)) +
    xlab(paste("log2([", conditions[1 + nb_rep[1]], "]/[", conditions[1], "])")) +
    ylab(expression(-log[10](" p-value "))) +
    scale_y_continuous(trans = "log1p")

  volcano_value <- volflu
  volcano_loc <- officer::ph_location(left = .5, top = 1.8, height = 4.5, width = 5)
  if (ibaq) {
    datib = ibaq_table(liste_prot_p1a2, datp1a2, aviBAQp1a2)
    datib2 = ibaq_table(liste_prot_a1p2, data1p2, aviBAQa1p2)
    xxp = paste0(nrow(datib), " proteins only quantified\n in ", conditions[1])
    xxp2 = paste0(nrow(datib2), " proteins only quantified\n in ", conditions[1 + nb_rep[1]])
    ylim_ib = c(min(c(datib$average_iBAQ, datib2$average_iBAQ)), max(c(datib$average_iBAQ, datib2$average_iBAQ)))
    if (!is.null(datib)) p2 = ibaq_barplot(datib, xxp, cols[2], ylim_ib)
    if (!is.null(datib2)) p3 = ibaq_barplot(datib2, xxp2, cols[3], ylim_ib)

    if ((!is.null(datib)) & (!is.null(datib2))) {
      volcano_value <- ggplotify::as.ggplot(gridExtra::arrangeGrob(p2, volflu, p3, nrow = 1, widths = c(1.5, 4, 1.5)))
    } else if ((!is.null(datib2)) & (is.null(datib))) {
      volcano_value <- ggplotify::as.ggplot(gridExtra::arrangeGrob(volflu, p3, nrow = 1, widths = c(4, 1.5)))
    } else if ((is.null(datib2)) & (!is.null(datib))) {
      volcano_value <- ggplotify::as.ggplot(gridExtra::arrangeGrob(p2, volflu, nrow = 1, widths = c(1.5, 4)))
    }
    if (!is.null(datib) || !is.null(datib2))
      volcano_loc <- officer::ph_location(left = .5, top = 1.8, height = 4, width = 8)
  }
  doc <- officer::ph_with(doc, value = volcano_value, location = volcano_loc)
  editable(volcano_value)

  # SAINT slides
  if (isTRUE(p$saint)) {
    fp_saint <- officer::fp_text(font.size = 16)
    parsSaint <- officer::block_list(officer::fpar(officer::ftext(paste0(
      "The shell command that has been used is:\n SAINTexpress-int ", "inter_", p$cond1_name, "_vs_", p$cond2_name, ".dat ",
      "prey_", p$cond1_name, "_vs_", p$cond2_name, ".dat ", "bait_", p$cond1_name, "_vs_", p$cond2_name, ".dat"), fp_saint)))
    parsinter <- officer::block_list(officer::fpar(officer::ftext(paste0(
      "Number of proteins with AvgP superior to ", p$saint_avgp, " and log2(fold-change) superior to ", log2(p$fold_change), ": ",
      as.character(length(saint_res$interestingData[, 4]))), fp_saint)))
    doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
    doc <- slide_title(doc, "SAINT results")
    doc <- officer::ph_with(doc, value = parsSaint, location = officer::ph_location(left = 0.5, top = 0.5, width = 10))
    doc <- officer::ph_with(doc, value = parsinter, location = officer::ph_location(left = 0.5, top = 1.5, width = 10))
    doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
    doc <- slide_title(doc, "SAINT results")
    doc <- officer::ph_with(doc, value = saint_res$plot, location = officer::ph_location(left = 0.5, top = 1.2, height = 6, width = 9))
    editable(saint_res$plot)
  }

  #################################
  # Result tables
  #################################
  step("Writing result tables")
  datt = data.frame(FCs, pvalu, res_abh$adjp[, 2], dat_filtered_FC)
  colnames(datt)[1:3] = c(paste("log2([ Mean in", conditions[1 + nb_rep[1]], "]/[ Mean in", conditions[1], "])"), "pvalue", "Adjusted pvalue")
  datt = datt[order(datt[, 3]), ]
  fc_col <- paste("log2([ Mean in", conditions[1 + nb_rep[1]], "]/[ Mean in", conditions[1], "])")

  # proteins with |log2FC| below the threshold (not tested for FDR)
  datFCinf = fill_missing_columns(dat_filtered_FC_init[which(abs(FC) < FCT), ])
  datFCinf_fin = data.frame(datFCinf$Fasta.headers, datFCinf$Majority.protein.IDs, datFCinf$Peptides, datFCinf$Mol..weight..kDa.,
                            FC[which(abs(FC) < FCT)], pval[which(abs(FC) < FCT)], rep(NA, length(pval[which(abs(FC) < FCT)])),
                            as.matrix(datFCinf[, grep(int_pattern, names(datdel), value = TRUE)]))
  colnames(datFCinf_fin)[1:7] = c("Fasta.headers", "Majority.protein.IDs", "Peptides", "Mol. weight kDa.", fc_col, "p.value", "Adjusted p.value")
  ann = annotate_proteins(datFCinf$Protein.IDs, datFCinf_fin$Fasta.headers)
  datFCinf_fin = add_annotation_columns(datFCinf_fin, ann)
  openxlsx::addWorksheet(wb, paste0("FC_inf_to", 2^FCT))
  openxlsx::writeDataTable(wb, paste0("FC_inf_to", 2^FCT), x = add_counts(datFCinf_fin, datFCinf$id))

  # proteins with |log2FC| above the threshold but not significant
  datnonsig = fill_missing_columns(datt[which(datt[, 3] > FDRT), ])
  datnonsig_fin = data.frame(datnonsig$Fasta.headers, datnonsig$Majority.protein.IDs, datnonsig$Peptides, datnonsig$Mol..weight..kDa.,
                             datt[which(datt[, 3] > FDRT), 1], datt[which(datt[, 3] > FDRT), 2], datt[which(datt[, 3] > FDRT), 3],
                             as.matrix(datnonsig[, grep(int_pattern, names(datdel), value = TRUE)]))
  colnames(datnonsig_fin)[1:7] = c("Fasta.headers", "Majority.protein.IDs", "Peptides", "Mol. weight kDa.", fc_col, "p.value", "Adjusted p.value")
  if (length(datnonsig_fin$Fasta.headers) != 0) {
    ann = annotate_proteins(datnonsig$Protein.IDs, datnonsig_fin$Fasta.headers)
    datnonsig_fin = add_annotation_columns(datnonsig_fin, ann)
  }
  openxlsx::addWorksheet(wb, paste0("NonSignif_with_FC_sup_to", 2^FCT))
  openxlsx::writeDataTable(wb, paste0("NonSignif_with_FC_sup_to", 2^FCT), x = add_counts(datnonsig_fin, datnonsig$id))

  # more abundant in condition 2: absent in 1 / present in 2 + significant with FC > 0
  p1a2 = dat[dat$id %in% liste_prot_a1p2, ]
  p1a2 = cbind(matrix(NA, nrow(p1a2), 3), p1a2)
  lFCp = datt[datt[, 1] > 0, ]
  lAD = lFCp[lFCp[, 3] < FDRT, ]
  colnames(p1a2) = colnames(lAD)
  liste_prot_12 = rbind(p1a2, lAD)
  np1a2 = nrow(p1a2)
  nlAD1 = nrow(lAD)
  mab12 <- more_abundant_table(liste_prot_12, int_pattern)
  datta <- add_counts(mab12$table, mab12$liste$id)

  # more abundant in condition 1
  a1p2 = dat[dat$id %in% liste_prot_p1a2, ]
  a1p2 = cbind(matrix(NA, nrow(a1p2), 3), a1p2)
  lFCp = datt[datt[, 1] < 0, ]
  lAD = lFCp[lFCp[, 3] < FDRT, ]
  colnames(a1p2) = colnames(lAD)
  liste_prot_21 = rbind(a1p2, lAD)
  na1p2 = nrow(a1p2)
  nlAD2 = nrow(lAD)
  mab21 <- more_abundant_table(liste_prot_21, int_pattern)
  datta2 <- add_counts(mab21$table, mab21$liste$id)

  # heatmaps of the proteins of interest (one slide)
  if (nrow(datta) != 0 || nrow(datta2) != 0) {
    doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
    doc = slide_title(doc, "The proteins of interest")
  }
  if (nrow(datta) != 0) {
    openxlsx::addWorksheet(wb, paste0("MoreAbIn", shortname2, "ThanIn", shortname1))
    openxlsx::writeDataTable(wb, paste0("MoreAbIn", shortname2, "ThanIn", shortname1), x = datta)
    ggp <- interest_heatmap(datta, c(vars1, vars2), mab12$liste$Gene.name, ldat2,
                            paste("More abund. in ", levels(conditions)[2], "\n than in", levels(conditions)[1]))
    doc <- officer::ph_with(doc, value = ggp$plot, location = officer::ph_location_left())
    editable(ggp$plot)
  }
  if (nrow(datta2) != 0) {
    openxlsx::addWorksheet(wb, paste0("MoreAbIn", shortname1, "ThanIn", shortname2))
    openxlsx::writeDataTable(wb, paste0("MoreAbIn", shortname1, "ThanIn", shortname2), x = datta2)
    ggp <- interest_heatmap(datta2, c(vars1, vars2), mab21$liste$Gene.name, ldat2,
                            paste("More abund. in ", levels(conditions)[1], "\n than in", levels(conditions)[2]))
    doc <- officer::ph_with(doc, value = ggp$plot, location = officer::ph_location_right())
    editable(ggp$plot)
  }

  # barplot of the number of proteins of interest
  doc <- officer::add_slide(doc, layout = "Two Content", master = "Office Theme")
  doc = slide_title(doc, "The proteins of interest")
  df_cumsum = data.frame(
    dd = as.factor(c(rep(paste("More abund. in", levels(conditions)[2], "\n than in", levels(conditions)[1]), 2),
                     rep(paste("More abund. in", levels(conditions)[1], "\n than in", levels(conditions)[2]), 2))),
    cat = as.factor(c("Present / Absent", "Differentially abundant", "Present / Absent", "Differentially abundant")),
    nb = c(np1a2, nlAD1, na1p2, nlAD2),
    label_ypos = c(np1a2 / 2, nlAD1 / 2 + np1a2, na1p2 / 2, nlAD2 / 2 + na1p2),
    label_xpos = c(1, 1.1, 2, 2.1))
  df_cumsum$dd = factor(df_cumsum$dd, levels = c(paste("More abund. in", levels(conditions)[2], "\n than in", levels(conditions)[1]),
                                                 paste("More abund. in", levels(conditions)[1], "\n than in", levels(conditions)[2])))
  df_cumsum$cat = factor(df_cumsum$cat, levels = c("Differentially abundant", "Present / Absent"))
  ggp = ggplot(data = df_cumsum, aes(x = dd, y = nb, fill = cat)) +
    geom_bar(stat = "identity") +
    geom_text(aes(x = label_xpos, y = label_ypos, label = nb), vjust = 1.6, color = "black", size = 3.5) +
    scale_fill_manual(values = c("#79a7fc", "#FFFFB3")) +
    xlab("") + ylab("Nb of proteins") +
    coord_flip() +
    theme(legend.title = element_blank(), legend.position = "bottom")
  doc <- officer::ph_with(doc, value = ggp, location = officer::ph_location(left = 1, top = 2, height = 5, width = 7))
  editable(ggp)

  #################################
  # Final files
  #################################
  step("Writing final files")
  openxlsx::addWorksheet(wb, "MoreAbundInOneOrTheOther")
  if (nrow(datta2) == 0) {
    openxlsx::writeDataTable(wb, "MoreAbundInOneOrTheOther", x = datta)
  } else if (nrow(datta) == 0) {
    openxlsx::writeDataTable(wb, "MoreAbundInOneOrTheOther", x = datta2)
  } else {
    openxlsx::writeDataTable(wb, "MoreAbundInOneOrTheOther", x = rbind(datta, datta2))
  }

  if (isTRUE(p$spectronaut)) {
    openxlsx::addWorksheet(wb, "SpectronautDataAfterConversion")
    openxlsx::writeDataTable(wb, "SpectronautDataAfterConversion", x = dat)
  }

  openxlsx::addWorksheet(wb, "DataVolcano")
  openxlsx::writeDataTable(wb, "DataVolcano", x = rbind(datFCinf_fin[, 1:10], datnonsig_fin[, 1:10], datta[, 1:10], datta2[, 1:10]))
  openxlsx::addWorksheet(wb, "Imputation")
  openxlsx::writeDataTable(wb, "Imputation", x = imput_matrix)
  openxlsx::addWorksheet(wb, "ParametersAnalysis")
  openxlsx::writeDataTable(wb, "ParametersAnalysis", x = parameters)
  openxlsx::addWorksheet(wb, "SampleCorrelation")
  openxlsx::writeDataTable(wb, "SampleCorrelation", x = data.frame(matrice_corr))
  openxlsx::addWorksheet(wb, "SamplCorrAfterNormImput")
  openxlsx::writeDataTable(wb, "SamplCorrAfterNormImput", x = data.frame(matrice_corri))

  if (isTRUE(p$conditional_pvalues)) {
    dcp = fill_missing_columns(dat_filtered_FC)
    ann = annotate_proteins(dcp$Protein.IDs, dcp$Fasta.headers)
    cond_tab = data.frame(dcp$Majority.protein.IDs, FCs, pvalu, cutoff[sel_fc], pval_cond[sel_fc], res_abh$adjp[, 2])
    colnames(cond_tab) = c("Majority.protein.IDs", fc_col, "p.value", "p.value at |log2FC| = threshold",
                           "Conditional p.value", "Adjusted p.value")
    cond_tab = add_annotation_columns(cond_tab, ann)
    cond_tab = cond_tab[order(cond_tab[, 9]), ]
    openxlsx::addWorksheet(wb, "ConditionalPvalues")
    openxlsx::writeDataTable(wb, "ConditionalPvalues", x = cond_tab)
  }

  add_description_sheet(wb, p, name1 = levels(conditions)[1], name2 = levels(conditions)[2],
                        fc_col = fc_col, FCT = FCT, FDRT = FDRT, ftxt = ftxt, count_cols = count_cols)

  xlsx_file <- out(paste0("Final_", levels(conditions)[2], "_vs_", levels(conditions)[1], ".xlsx"))
  openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)
  pptx_file <- out(paste0("Differential_Analysis_", levels(conditions)[2], "_vs_", levels(conditions)[1], ".pptx"))
  print(doc, target = pptx_file)

  # parameters and R session, to re-run the analysis
  p_save <- p
  p_save$input_file <- if (is.null(p$input_file)) ds$chemin else p$input_file
  p_save$delete_columns <- del_cols
  param_file <- out("parameters.yml")
  write_params(p_save, param_file)
  session_file <- out("sessionInfo.txt")
  writeLines(c(paste("pairwiseApp", pairwise_version()),
               paste("Input file:", ds$chemin, "(", ds$sizee, "bytes ), md5:",
                     if (!is.null(p$input_file) && file.exists(p$input_file)) unname(tools::md5sum(p$input_file)) else "NA"),
               "", utils::capture.output(utils::sessionInfo())), session_file)

  files <- c(xlsx_file, pptx_file, param_file, session_file, overlap_files, saint_res$files)
  if (isTRUE(p$editable_pptx) && file.exists(editable_file)) files <- c(files, editable_file)

  Tf = proc.time()
  log(paste0("Time to execute : ", round((Tf - Ts)[3], 2), " secs"))

  invisible(list(
    files = files,
    counts = c(initial = n0, after_deletion = n3, after_filter1 = n4, after_filter2 = n5,
               present_absent_cond1 = length(liste_prot_p1a2), present_absent_cond2 = length(liste_prot_a1p2),
               significant_up = nlAD1, significant_down = nlAD2),
    results = datt, imputed = imput_matrix, overlap = overlap_results, saint = saint_res))
}


#################################
# Helpers used by run_pairwise() (plot code moved out of the main function)
#################################

pairwise_version <- function() "1.0.0"

# p-value that each protein would have with |log2FC| equal to the threshold, given
# its moderated standard error (limma fit, last coefficient = tested contrast)
conditional_cutoff <- function(fit, threshold) {
  j = ncol(fit$p.value)
  2 * pt(-threshold / (fit$stdev.unscaled[, j] * sqrt(fit$s2.post)), df = fit$df.total)
}

# "Description" sheet (first sheet of the workbook): content of each sheet and
# meaning of the main columns
add_description_sheet <- function(wb, p, name1, name2, fc_col, FCT, FDRT, ftxt, count_cols) {
  s1 <- short_name(name1); s2 <- short_name(name2)
  fc_txt <- paste0("|log2(fold-change)| > ", round(FCT, 4), " (fold-change ", p$fold_change, ")")
  sheets <- c(
    "Description" = "This sheet: content of each sheet of the workbook and meaning of the main columns.",
    "without_enough_values" = paste0("Proteins removed: neither tested nor present/absent. ", ftxt$definition, " ", ftxt$tested, " ", ftxt$present,
                                     " The Reason column gives the reason for each protein; the intensities are those of the input file."),
    "SAINT" = paste0("SAINTexpress output (one row per prey): Bait = name of condition 2; Prey = protein identifier; PreyGene = gene name (from the input file); ",
                     "Intensity / ctrlIntensity = intensities in the bait / control samples (separated by |); AvgP = average probability of interaction over the bait replicates; ",
                     "MaxP = maximal probability; SaintScore; OddsScore; FoldChange = SAINTexpress fold-change between bait and control intensities; ",
                     "BFDR = Bayesian false discovery rate."),
    "Bait" = "Bait file given to SAINTexpress: sample name, bait name (name of condition 2, or Ctrl for the controls), T (test) or C (control).",
    setNames(paste0("Proteins present in ", name1, " and absent from ", name2, ". ", ftxt$present,
                    if (isTRUE(p$ibaq)) paste0(" 'Average iBAQ ", name1, "' = mean log2(iBAQ) over the samples of ", name1, ".") else ""),
             paste0("PresIn", s1, "_AbsIn", s2)),
    setNames(paste0("Proteins present in ", name2, " and absent from ", name1, ". ", ftxt$present,
                    if (isTRUE(p$ibaq)) paste0(" 'Average iBAQ ", name2, "' = mean log2(iBAQ) over the samples of ", name2, ".") else ""),
             paste0("PresIn", s2, "_AbsIn", s1)),
    setNames(paste0("Tested proteins whose |log2(fold-change)| is below the threshold (", round(FCT, 4), "): not considered for the FDR control (adjusted p-value = NA). ",
                    ftxt$tested),
             paste0("FC_inf_to", 2^FCT)),
    setNames(paste0("Proteins with ", fc_txt, " whose adjusted p-value is above the FDR level (", FDRT, ")."),
             paste0("NonSignif_with_FC_sup_to", 2^FCT)),
    setNames(paste0("Proteins more abundant in ", name2, " than in ", name1, ": present in ", name2, " and absent in ", name1,
                    " (no p-value), and significant proteins (adjusted p-value < ", FDRT, ", ", fc_txt, ") with a positive log2(fold-change)."),
             paste0("MoreAbIn", s2, "ThanIn", s1)),
    setNames(paste0("Proteins more abundant in ", name1, " than in ", name2, ": present in ", name1, " and absent in ", name2,
                    " (no p-value), and significant proteins (adjusted p-value < ", FDRT, ", ", fc_txt, ") with a negative log2(fold-change)."),
             paste0("MoreAbIn", s1, "ThanIn", s2)),
    "MoreAbundInOneOrTheOther" = "The two previous sheets together.",
    "SpectronautDataAfterConversion" = "Spectronaut input data after conversion to the column names used by the application (Intensity.<sample>, Peptides.<sample>, ...).",
    "DataVolcano" = "First columns (identifiers, log2(fold-change), p-value, adjusted p-value) of all the proteins of the sheets FC_inf_to, NonSignif_with_FC_sup_to and MoreAbIn: data of the volcano plot.",
    "Imputation" = paste0("log2 intensities after filter 1, normalisation (", p$normalization, ") and imputation (", p$imputation,
                          "): values used for the statistical test of the tested proteins. The sheet also contains the present/absent proteins, which are not tested (they can keep missing values in the condition where they are absent)."),
    "ParametersAnalysis" = "Parameters of the analysis.",
    "SampleCorrelation" = "Pearson correlation between samples computed on log2 intensities (pairwise complete observations), before filter 1, normalisation and imputation.",
    "SamplCorrAfterNormImput" = "Pearson correlation between samples computed on log2 intensities after filter 1, normalisation and imputation.",
    "ConditionalPvalues" = paste0("Proteins with ", fc_txt, ": p-value, p-value that the protein would have with |log2(fold-change)| equal to the threshold, ",
                                  "conditional p-value (ratio of the two) and adjusted p-value computed from the conditional p-values.")
  )
  present <- openxlsx::sheets(wb)
  sheets <- sheets[names(sheets) %in% c("Description", present)]
  sheets <- sheets[order(match(names(sheets), c("Description", present)))]

  columns <- c(
    "Uniprot accession / Protein name / Gene name" = "Extracted from the protein identifiers and FASTA headers of the input file.",
    "Majority.protein.IDs, Fasta.headers, Peptides, Mol..weight..kDa." = "Columns of the input file (Spectronaut: gene names are in Majority.protein.IDs).",
    setNames(paste0("Difference between the mean log2 intensities of ", name2, " and ", name1, " (after normalisation and imputation)."), fc_col),
    "pvalue / p.value" = paste0("p-value of the ", if (p$test == "LIMMA") "moderated t-test (limma)" else "paired moderated t-test (limma)", "."),
    "Adjusted pvalue / Adjusted p.value" = paste0("p-value adjusted with the adaptive Benjamini-Hochberg procedure (proportion of non-differential proteins estimated with the ",
                                                  p$pi0_method, " method of cp4p) among the proteins with ", fc_txt,
                                                  if (isTRUE(p$conditional_pvalues)) ", computed from the conditional p-values (see ConditionalPvalues)" else "",
                                                  ". A protein is significant if its adjusted p-value is below ", FDRT, "."),
    "Intensity.<sample> / LFQ.intensity.<sample>" = "Intensities of the input file (not log-transformed, not normalised), except in the Imputation sheet (log2, normalised, imputed).",
    "iBAQ.<sample>" = "iBAQ values of the input file.",
    setNames(paste0("Number of quantified values of the protein in each condition (", ftxt$tested, ")"),
             paste(count_cols, collapse = " / ")),
    "Reason (sheet without_enough_values)" = "Why the protein is neither tested nor present/absent."
  )

  openxlsx::addWorksheet(wb, "Description")
  openxlsx::writeData(wb, "Description", "Content of the sheets", startRow = 1, startCol = 1)
  openxlsx::writeDataTable(wb, "Description", data.frame(Sheet = names(sheets), Content = unname(sheets)),
                           startRow = 2, tableName = "Sheets")
  r <- length(sheets) + 5
  openxlsx::writeData(wb, "Description", "Main columns", startRow = r - 1, startCol = 1)
  openxlsx::writeDataTable(wb, "Description", data.frame(Column = names(columns), Meaning = unname(columns)),
                           startRow = r, tableName = "Columns")
  openxlsx::setColWidths(wb, "Description", cols = 1:2, widths = c(45, 150))
  openxlsx::addStyle(wb, "Description", openxlsx::createStyle(wrapText = TRUE, valign = "top"),
                     rows = 2:(r + length(columns)), cols = 1:2, gridExpand = TRUE)
  openxlsx::addStyle(wb, "Description", openxlsx::createStyle(textDecoration = "bold", fontSize = 12),
                     rows = c(1, r - 1), cols = 1)
  n <- length(openxlsx::sheets(wb))
  openxlsx::worksheetOrder(wb) <- c(n, seq_len(n - 1))
  invisible(wb)
}

# Heatmap of a correlation (or Jaccard) matrix between samples
corr_heatmap <- function(z, gradientRate) {
  text <- element_text(colour = "black", size = 8, face = "bold")
  suppressWarnings(ggplot(reshape2::melt(z), aes(x = Var1, y = Var2, fill = value)) + geom_tile() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), legend.text = text) +
    theme(legend.title = element_blank()) + labs(x = "", y = "") +
    scale_fill_gradientn(colours = colorRampPalette(c("white", "lightblue", "darkblue"))(101),
                         values = c(pexp(seq(0, 1, 0.01), rate = gradientRate), 1), limits = c(0, 1)))
}

# Dendrogram of samples coloured by condition
sample_dendrogram <- function(hv, coln, conditions, title, colour = TRUE) {
  ddhv <- as.dendrogram(hv)
  lddhv <- labels(ddhv)
  condhc = rep(0, length(lddhv))
  for (lb in 1:length(lddhv)) {
    for (lc in 1:length(coln)) {
      if (lddhv[lb] == coln[lc]) { condhc[lb] = conditions[lc] }
    }
  }
  colv = as.numeric(as.factor(condhc))
  colv[colv == 1] = "coral1"
  colv[colv == 2] = "cyan4"
  if (!colour) colv = "black"
  ggd = ggdendro::ggdendrogram(hv, rotate = TRUE, theme_dendro = FALSE)
  ggd = suppressWarnings(ggd + theme(axis.text.y = element_text(hjust = 1, colour = colv)))
  ggd + xlab("Sample") + ylab("Height") + ggtitle(title)
}

# Small text block at the bottom left of a slide
note_block <- function(doc, line1, line2) {
  fp_normal <- officer::fp_text(font.size = 12)
  pars <- officer::block_list(officer::fpar(officer::ftext(line1, fp_normal)),
                              officer::fpar(officer::ftext(" ", fp_normal)),
                              officer::fpar(officer::ftext(line2, fp_normal)))
  officer::ph_with(x = doc, value = pars, location = officer::ph_location(left = 0.2, top = 6.2, height = 1.2, width = 4))
}

# Jaccard index between samples computed on the missing-value pattern
jaccard_missing <- function(df, margin) {
  if (margin == 1 | margin == 2) {
    res = matrix(0, ncol(df), ncol(df))
    for (i in 1:ncol(df)) {
      for (j in 1:ncol(df)) {
        M_00 <- apply(cbind(df[, i], df[, j]), margin, sum) == 0
        M_11 <- apply(cbind(df[, i], df[, j]), margin, sum) == 2
        if (margin == 1) {
          dff <- cbind(df[, i], df[, j])[!M_00, ]
          JSim <- sum(M_11) / nrow(dff)
        } else {
          dff <- cbind(df[, i], df[, j])[, !M_00]
          JSim <- sum(M_11) / length(dff)
        }
        res[i, j] = JSim
      }
    }
    colnames(res) = colnames(df)
    rownames(res) = colnames(res)
    return(res)
  }
}

# Mean log2(iBAQ) over the samples of a condition (present/absent proteins)
mean_log2_ibaq <- function(tab, sample_names) {
  m = as.matrix(tab[, grep("iBAQ", names(tab), value = TRUE)])
  ibaqnames = paste("iBAQ.", sample_names, sep = "")
  if (nrow(m) > 1) {
    if (!is.numeric(m[, 1])) m = apply(gsub(" ", "", m), 2, function(x) as.numeric(x))
    m = log2(as.matrix(m[, ibaqnames]))
  } else {
    if (!is.numeric(m[, 1])) m[] = as.numeric(gsub(" ", "", m))
    mode(m) = "numeric"
    m = log2(t(as.matrix(m[, ibaqnames])))
  }
  m[is.infinite(m)] = NA
  apply(m, 1, mean, na.rm = T)
}

ibaq_table <- function(ids, tab, avibaq) {
  if (length(ids) > 1) {
    d = data.frame(names = tab[, 1], average_iBAQ = as.numeric(avibaq))
    d = d[which(!is.na(d[, 2])), ]
    d$average_iBAQ = as.numeric(d$average_iBAQ)
    d$names = as.factor(d$names)
    d$names = factor(d$names, levels = d$names[order(d$average_iBAQ)])
    d
  } else if (length(ids) == 1) {
    data.frame(names = tab$`Uniprot accession`[1], average_iBAQ = mean(avibaq, na.rm = T))
  } else NULL
}

ibaq_barplot <- function(d, xlabel, colour, ylim_ib) {
  ggplot(d, aes(x = names, y = average_iBAQ)) +
    ylab(expression(log[2](" iBAQ "))) + xlab(xlabel) +
    geom_bar(stat = "identity", fill = colour) +
    ggtitle(" ") + theme_classic() +
    coord_flip(ylim = ylim_ib) +
    theme(axis.text.x = element_text(face = "bold", colour = "black", size = 10)) +
    theme(axis.text.y = element_blank()) +
    theme(axis.title = element_text(size = 15, face = "bold"))
}

# Table of proteins more abundant in a condition (present/absent + significant)
more_abundant_table <- function(liste, int_pattern) {
  code = liste$Fasta.headers
  coln = as.character(code)
  if (length(liste$Fasta.headers) != 0) {
    for (i in 1:length(coln)) {
      car = unlist(strsplit(coln[i], "[|]"))
      coln[i] = car[2]
    }
    liste = cbind(liste, coln)
    colnames(liste)[ncol(liste)] = "Gene.name"
  }
  liste = fill_missing_columns(liste)
  datta = data.frame(liste$Fasta.headers, liste$Majority.protein.IDs, liste$Peptides, liste$Mol..weight..kDa.,
                     as.matrix(liste[, grep("log2", names(liste), value = TRUE)]), liste$pvalue, liste$`Adjusted pvalue`,
                     as.matrix(liste[, grep(int_pattern, names(liste), value = TRUE)]))
  if (nrow(datta) != 0) {
    colnames(datta)[c(1:7)] = c("Fasta.headers", "Majority.protein.IDs", "Peptides", "Mol. weight kDa.", colnames(liste)[1], "p.value", "Adjusted p.value")
    ann = annotate_proteins(liste$Protein.IDs, datta$Fasta.headers)
    datta = add_annotation_columns(datta, ann)
  }
  list(table = datta, liste = liste)
}

# Heatmap of log2 intensities of the proteins of interest
interest_heatmap <- function(datta, sample_cols, gene_names, ldat2, title) {
  distCor <- function(x) as.dist(1 - cor(x, use = "pairwise"))
  zClust <- function(z, method = "complete") {
    a = distCor(t(z)); a[is.na(a)] = 0
    hcl_row <- hclust(a, method = method)
    list(Clust.row = hcl_row)
  }
  da = as.matrix(datta[, sample_cols])
  for (i in 1:ncol(da)) { da[, i] = gsub(",", ".", da[, i]) }
  da = apply(da, 2, as.numeric)
  da = log2(da)
  da[da < 0] = 0
  if (is.vector(da)) { da = t(as.data.frame(da)) }
  rownames(da) = gene_names

  clustered <- nrow(da) > 3
  pd <- as.data.frame(cbind(rownames(da), da))
  colnames(pd)[1] = "Proteins"
  pd.m = reshape2::melt(pd, id.vars = "Proteins", variable.name = "Samples")
  if (clustered) {
    variance_colonnes <- apply(da, 2, var)
    ord <- zClust(da[, variance_colonnes > 0], method = "ward.D")$Clust.row$order
    for (i in seq_along(pd.m$Proteins)) {
      if (grepl(";", pd.m$Proteins[i])) pd.m$Proteins[i] <- sub(";.*", "", pd.m$Proteins[i])
    }
    pd.m$value = as.numeric(pd.m$value)
    pd.m$Proteins <- factor(pd.m$Proteins, levels = make.unique(pd.m$Proteins[ord]))
  } else {
    pd.m$value = as.numeric(pd.m$value)
    pd.m$Proteins <- as.factor(pd.m$Proteins)
  }
  pd.m = as.data.frame(pd.m)
  plot = ggplot(pd.m, aes(x = Samples, y = Proteins, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colours = c("grey80", "yellow", "red"), values = c(0, 0.4, 1), limits = c(0, max(ldat2[, 2:ncol(ldat2)], na.rm = T))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size = 2)) +
    labs(fill = "log2(int.)") +
    ggtitle(title)
  list(plot = plot)
}
