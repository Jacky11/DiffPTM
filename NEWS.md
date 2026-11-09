# DiffPTM 2.0.0 (development)

* Default change: the usage (DPU) model now pairs a site with its protein through a fixed sample
  effect when the two assays share sample names (`blocking = "auto"` is the *per-site paired single
  model*, equivalent to limma on the within-sample site-minus-protein differences). The consensus
  model (random sample effect, `blocking = "random"`) pools the residual variances of the site and of
  its protein and shares one correlation between sites; in the real datasets of the benchmark the
  residual variance of a site was 3-9 times its protein's and this ratio varied between sites, and in
  simulations with such variances the consensus model was anti-conservative while the per-site model
  stayed closest to the nominal level. The consensus model remains available.
* New `dpu_diagnostics()`: within-sample correlation, site/protein variance ratio and its spread
  over sites, replicates and prior degrees of freedom, with a recommended pairing.
  The spread is the excess of the observed variance of the log variance ratio over its sampling
  expectation, which accounts for the residual correlation of site and protein (Monte Carlo with a
  fixed internal seed; the user's random stream is unchanged); `rho_obs_pooled` and
  `log_var_ratio_excess` are reported.
* New `usage_beta_check()`: refits the usage model with the protein term scaled by `beta`
  (default 0.8 and 1.2) and flags the calls that are not robust to a departure from the unit
  protein coefficient; also returns the protein log2 fold change of each site.
* `method_limma(blocking = "random")` stops with an explicit message when the consensus
  correlation cannot be estimated (too few residual degrees of freedom).
* `simulate_ptm_features(site_sd_factor = )`: site sample effects with a larger SD than the
  protein's (default 1, unchanged).
* `methods_paragraph()` and `write_results_xlsx()` describe the pairing that was used (`blocking`).
* Default change: the consensus-correlation usage model (`method_limma()`, sample as a random
  effect) is now fitted on complete site-protein pairs only (`complete_pairs = TRUE`). With the
  previous all-observation fit, a sample with only the site or only the protein informed the
  other component through the shared correlation; in simulations this made the test
  anti-conservative for sites with unpaired values (mean FDP up to 0.12 at a nominal 0.05),
  whereas the complete-pairs fit was calibrated like the per-site paired model
  (`blocking = "fixed"`). `complete_pairs = FALSE` restores the previous behaviour.
* Default change: `normalise()` now centres each assay on the global median of all samples
  (`"global_median"`, as MSstats and MSstatsTMT). The previous default, per-condition median
  centering (`"median"`), preserves a difference of mean loading between conditions, which the
  usage test attributes to every site; it remains available.
* `simulate_ptm_features()` gains options for sensitivity analyses: site clustering on proteins,
  condition-specific correlation, heterogeneous variances, unmodelled batches, PTM loading offsets,
  outlying values, unbalanced designs (`n_rep = c(n_A, n_B)`) and unbalanced missingness; the
  default simulation is unchanged.

* Result tables now carry the degrees of freedom (`df`) of each method's test and a 95 % confidence
  interval (`ci_lo`, `ci_hi`; `add_confidence_interval()`), exported in the Excel workbook and shown in
  the app; `method_msqrob2ptm(robust = FALSE)` fits ordinary least squares instead of M-estimation.
* Naming: the single model with a random sample term (one consensus within-sample correlation)
  is the *consensus single model*; `blocking = "fixed"` gives the *per-site paired single model*;
  `pair_by_sample = FALSE` the *unpaired single model*. The article and the app use these names.
* Documentation of the single model corrected: when the two assays share sample names,
  the sample term is fitted as a *random* effect with one consensus within-sample correlation
  (`limma::duplicateCorrelation`), not as a fixed effect; the fixed effect (exact per-site pairing,
  as in msqrob2PTM) is `blocking = "fixed"`. `method_limma(moderate = FALSE)` returns ordinary
  t-statistics (used in the benchmark to separate pairing from moderation).
* `simulate_ptm_features()` accepts a per-site `rho` and a protein coefficient `beta` (sites
  following their protein with a coefficient other than one); new benchmark scripts
  `inst/benchmark/run_factorial.R` (pairing x moderation with one estimator) and
  `run_sensitivity.R` (unit coefficient, heterogeneous correlation, summarisation filter) and
  `run_combination.R` (per-site calls of all strategies, to evaluate intersection/union/majority rules),
  `run_extra_seeds.R` (additional MSstatsPTM seeds; robust vs non-robust msqrob2PTM).
* Importers checked on real MaxQuant (combined search) and Spectronaut (pivot report) exports and
  on example exports of Proteome Discoverer, FragPipe (label-free and TMT) and DIA-NN
  (`inst/extdata/formats/`, built by `inst/scripts/make_format_fixtures.R`); one vignette per
  software plus an overview of the interface.
* Fixed: MaxQuant intensities above 2^31 were read as `integer64` and transformed to nonsense;
  Spectronaut exports with semicolons and decimal commas were read as missing.
* `import_table()` detects the decimal comma; `check_scale()` warns when the scale of the values
  does not match `log2_transform`; `rename_samples()` / `keep_samples()` select and rename the
  samples of an assay (combined searches, re-acquired runs); `guess_annotation()` accepts token
  ranges (`"1-2"`) and `-1` (all tokens but the last); `make_contrasts(type = "vs_reference")`.
* `locate_sites()` fills residues and protein positions of peptidoform tables (DIA-NN precursors,
  FragPipe modified peptides) from the modified sequence and the FASTA; Proteome Discoverer sites
  are parsed from "Modifications in Master Proteins"; Spectronaut `PTM.SiteAA`, `PTM.Multiplicity`
  and `PTM.FlankingRegion` are read; a short `label` (gene + site) is built for every feature.
* App: new title; design step with sample-name patterns, token ranges and exclusion of unassigned
  samples; example dataset shipped (`inst/extdata/example`); processing step merges redundant
  site rows and locates sites of peptide tables; QC adds features observed per sample and intensity
  distributions; reference-based contrasts; results page with gene/site labels, counts per method
  and a filtered discordance table; a common plot theme (`theme_diffptm()`), PCA, counts and
  distribution plots (`plot_pca()`, `plot_counts()`, `plot_intensity_distributions()`,
  `plot_detected_per_sample()`).
* `filter_min_observed()`: `min_n` is now the second argument (`assay` moved after `min_pct`).

* The statistical approaches are now called "methods" everywhere (previously "engines"):
  `run_methods()`, `method_limma()`, `method_msstatsptm()`, `method_msstatsptm_native()`,
  `method_msqrob2ptm()`, `available_methods()`, `compare_methods()`; the result column is `method`.
* `merge_redundant_ptm()` merges site rows reported once per protein accession with identical
  quantities (Spectronaut PTM site reports).

* Rewritten as an R package; the Shiny app is a thin layer over the exported API.
* One data model (`ptm_experiment`) with per-assay sample annotation, `subject` for paired /
  blocked designs and a processing log that generates the Materials & Methods paragraph.
* Importers for MaxQuant, Spectronaut, Proteome Discoverer, FragPipe, DIA-NN and a generic
  wide format; column specifications are data and can be overridden in the app.
* PTM-to-protein matching by identifier or by peptide sequence (FASTA, I/L-equivalent),
  shared peptides averaged; no NA-to-zero conversion.
* Site collapsing and FragPipe-style multisite features.
* Three statistical methods with a common result schema: limma single model (DPA + DPU),
  MSstatsPTM, msqrob2PTM; method concordance and consensus calls.
* Missing values tolerated by the limma method (imputation optional); detection-profile
  discordance between PTM features and their parent protein.
* Reports: Excel with README + Methods, editable PowerPoint, PDF, session `.rds`.
