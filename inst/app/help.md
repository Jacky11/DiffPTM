# DiffPTM — help

## What the app does

A modified peptide's intensity changes when the modification changes **and** when its protein
changes. DiffPTM separates the two:

* **Differential PTM abundance (DPA)** — does the site (or modified peptide) change between
  conditions, whatever the reason?
* **Differential PTM usage (DPU)** — does the site change *relative to its parent protein*
  (is the stoichiometry of the modification altered)? This is the test that needs the
  total-proteome table.

Three statistical methods answer the DPU question and can be run side by side on the same
processed data:

| Method | How | Reference |
|---|---|---|
| Single model (limma) | site and protein intensities fitted together; the condition × component interaction is the usage change. Paired (uses the within-sample site–protein correlation) when PTM and proteome samples share the same names | Wolski et al., Methods Mol Biol 2027; Ritchie et al. 2015 |
| MSstatsPTM | site and protein modelled separately, then the protein fold change is subtracted and the uncertainties added | Kohler et al., MCP 2023 |
| msqrob2PTM | protein intensity subtracted from the site in each sample, then the differences are modelled (robust regression with variance moderation) | Demeulemeester et al., MCP 2024 |

Estimates of the usage change are the same for the three; their standard errors differ (pairing,
variance moderation), hence different numbers of significant sites. With few replicates and
paired samples the single model is the most sensitive; MSstatsPTM is the choice when the PTM
and proteome samples are not the same.

## The six steps

1. **Import** — select the software and upload the proteome and PTM tables (see formats
   below). Tick off "Intensities are raw" for tables already on the log2 scale (FragPipe
   TMT-Integrator); the app warns when the scale looks wrong.
2. **Design** — give every sample a condition. When the PTM and proteome samples of the
   same biological sample carry different names (`Phospho_mock_1` and `mock_1`, or
   `Phospho-01_NS_R1` and `Proteome-01_NS_R1`), enter a regular expression per table whose
   capture group is the common name; samples that do not match are dropped. "Pre-fill from
   sample names" builds conditions from the name tokens (`-1` = everything but the replicate
   number). Optional: a second factor (time, batch), a subject (paired samples) and the TMT
   mixture. Samples left without a condition are excluded.
3. **Processing** — remove reverse hits / contaminants, merge site rows repeated for
   several protein accessions (Spectronaut), filter on localization probability and on
   observed values, normalise, match sites to proteins (by identifier, or by locating the
   peptide in the FASTA, which also gives site positions for peptide tables), optionally
   collapse peptidoforms to sites. Imputation is optional and **off by default**: the
   methods fit on observed values, and imputing before the correction inflates false
   discoveries. The quality-control plots show the number of features observed per sample
   (failed runs stand out), the Jaccard index of detection between samples, a PCA and the
   intensity distributions.
4. **Statistics** — choose the methods, the analyses, the contrasts (each condition against a
   reference, every pair, or each condition against the others), the FDR and an optional
   fold-change threshold.
5. **Results** — volcano plot with gene/site labels, p-value histogram, PTM-vs-protein fold
   changes (points off the diagonal change in usage), counts per contrast and method,
   agreement between methods, the table of significant features, and the PTM features
   detected in samples where their protein is not (detection discordance).
6. **Download** — Excel workbook (all tables, a README sheet and the automatically written
   Materials & Methods paragraph), figures as editable PowerPoint and PDF, and the session as
   an `.rds` file that `readRDS()` reloads in R.

## Input formats

| Software | Proteome table | PTM table | Notes |
|---|---|---|---|
| MaxQuant | `proteinGroups.txt` (LFQ intensities when present) | `Phospho (STY)Sites.txt`, `GlyGly (K)Sites.txt`, ... | combined searches list the proteome runs in the site table too: rename the `Phospho_` runs in step 2; multiplicity columns `___1/2/3` are expanded, kept or summed |
| Spectronaut | protein pivot report (`PG.Quantity`) | PTM site pivot report (`PTM.Quantity`, `PTM.SiteProbability`, `PTM.SiteAA`, `PTM.SiteLocation`) | semicolon separators and decimal commas are detected; `Filtered` = missing; one row per accession of a protein group is merged in step 3 |
| Proteome Discoverer | Proteins export | PeptideGroups export | tab-separated text exports; per-file abundances (`Abundance: F1: Sample`), normalised ones preferred; sites read from "Modifications in Master Proteins" |
| FragPipe | `combined_protein.tsv` (label-free) or `abundance_protein_*.tsv` (TMT) | `combined_modified_peptide.tsv` (label-free) or `abundance_single-site_*.tsv` (TMT) | TMT-Integrator tables are log2 ratios: untick "Intensities are raw" |
| DIA-NN | `report.pg_matrix.tsv` | `report.pr_matrix.tsv` (precursors; collapse to sites in step 3 with the FASTA) or a site matrix | columns are named after the raw files; the `.raw` / `.d` extension is removed |
| Generic | wide table: `protein_id`, `gene`, `intensity_<sample>` | wide table: `feature_id`, `protein_ids`, `gene`, `stripped_sequence`, `modified_sequence`, `residue`, `position`, `localization_prob`, `ptm_type`, `intensity_<sample>` | any other software, after reshaping |

"Advanced: column mapping" accepts a regular expression with one capture group naming the
sample, for exports whose column names differ from the defaults. The package vignettes
(`browseVignettes("DiffPTM")`) walk through one dataset per software.

## Limitations

* MSstatsPTM and msqrob2PTM are optional dependencies; when they are not installed the
  corresponding options are disabled. MSstatsPTM takes minutes for thousands of sites.
* Sequence-based matching and site positions of peptide tables require the FASTA used for the
  searches.
* Multisite features are grouped by the set of sites only (no probability weighting).
