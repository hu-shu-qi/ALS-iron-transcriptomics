# Iron-related transcription in ALS spinal cord

Analysis code for *Cell-type and regional patterns of iron-related transcription in amyotrophic lateral sclerosis spinal cord*.

Version 1.1.0. Original code is provided under the MIT license.

## Contents

| Path | Purpose |
|---|---|
| `scripts/` | Numbered analysis entry points |
| `scripts/lib/` | Shared scoring, count-matrix and figure functions |
| `config/` | Public sample mappings, fixed nucleus annotations and feature definitions |
| `resources/iron_gene_sets.gmt` | Fourteen Human MSigDB v2026.1.Hs gene sets |
| `source_data/` | Aggregate estimates and fixed coordinates for Figures 2–6 and S1–S5 |
| `renv.lock` | R package versions |
| `LICENSE` | MIT terms for original code |

The analysis starts from public processed matrices and the study's fixed single-nucleus annotation. Clinical model fitting requires authorized access to the institutional workbook. The aggregate figure inputs can be used without patient-level records.

## Environment

Use R 4.5.1 and the package versions in `renv.lock`. The environment was tested on macOS arm64. R itself, compiler/system libraries and the Arial font must be available separately. The figures use Cairo, svglite and ragg for PDF, SVG, PNG and TIFF output.

Install renv 1.3.1 in an accessible bootstrap library. Run the following from the repository root, using new directories for the restored library and restoration workspace:

```sh
Rscript --vanilla scripts/1_restore_environment.R renv.lock /absolute/path/to/restore_work /absolute/path/to/R_library
```

Set the library and input/output locations for subsequent commands:

```sh
export R_LIBS="/absolute/path/to/R_library"
export R_LIBS_USER="$R_LIBS"
export R_LIBS_SITE="$R_LIBS"
export LC_ALL=en_US.UTF-8
ALS_INPUT="/absolute/path/to/public_inputs"
ALS_OUT="/absolute/path/to/results"
mkdir -p "$ALS_OUT"
```

Each analysis output directory supplied below must be new. Keep the generated RDS objects with their associated BPCells directories.

## Entry points

| Script | Analysis and principal output |
|---|---|
| `1_restore_environment.R` | Restore packages from `renv.lock` |
| `2_prepare_bulk.R` | Public normalized matrix → filtered log2 expression and donor metadata |
| `3_prepare_single_nucleus.R` | Public H5 counts + fixed annotation → annotated Seurat object and BPCells matrices |
| `4_aggregate_single_nucleus.R` | Counts → donor–cell-type pseudobulk profiles |
| `5_aggregate_spatial.R` | Spatial objects → donor–region count matrices |
| `6_differential_expression.R` | Bulk limma and cell-type edgeR models, candidates and cell-type support |
| `7_single_nucleus_scores.R` | Nucleus-level UCell, AUCell and singscore → donor–cell-type scores |
| `8_gene_set_models.R` | Bulk/spatial scoring and donor-level gene-set models across all three layers |
| `9_spatial_expression.R` | Regional edgeR models and five-model support |
| `10_spatial_lodo.R` | Leave-one-donor-out regional expression models |
| `11_spatial_qc_counts.R` | Reaggregate spatial counts under four spot-quality conditions |
| `12_spatial_qc_models.R` | Regional expression sensitivity to spot-quality filtering |
| `13_conditional_enrichment.R` | Two-stage gene-set enrichment, exact intervals and joint BH correction |
| `14_clinical_associations.R` | Spearman correlations, patient bootstrap intervals, HC3 models and diagnostics |
| `15_prepare_plot_data.R` | Explicit model results → figure input files |
| `16_render_figures.R` | Export Figures 2–6 and S1–S5 |
| `17_render_spatial_lodo.R` | Export Figure S3 from LODO results or its published aggregate summary |

Scripts in `scripts/lib/` are loaded by the entry points and remain part of the runnable package. Their individual purposes are:

| Module | Purpose |
|---|---|
| `bpcells.R` | Locate the backing count matrices when an annotated object has moved |
| `scoring.R` | Three scoring algorithms, linear models and multiple-testing correction |
| `plot_style.R` | Shared palettes, labels, plotting theme and figure exporters |
| `figure2.R`–`figure6.R` | Construct each main analysis figure |
| `figureS1.R` | Microglial three-algorithm estimates and confidence intervals |
| `supplementary_figures.R` | Figures S2, S4 and S5 |

## Public data and fixed resources

### Bulk RNA sequencing

Download [GSE287256 normalized counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287256/suppl/GSE287256_ALS_bulkRNAseq_results_bulk_counts_normalized.csv.gz).

`config/GSE287256_samples.csv` maps the twelve expression columns to GEO accessions, donors, group and sex. The input requires `gene_ID`, `gene_name` and the twelve `expression_sample_id` columns. Values are the original investigators' DESeq2-normalized counts. Duplicate symbols are selected by highest median expression, then gene identifier. Genes require expression >1 in at least three samples and nonzero variance. The analysis applies log2(x+1).

### Single-nucleus RNA sequencing

Download and extract [GSE287257_RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287257/suppl/GSE287257_RAW.tar). Put the twelve H5 matrices in one directory. `config/sn_raw_files.csv` gives their filenames, sizes and MD5 hashes.

The reconstruction uses:

- `sn_frozen_annotation.csv.gz`: 100,265 retained nuclei, cell types, group, QC fields, clustering labels and UMAP coordinates.
- `sn_frozen_features.csv`: the 32,869 retained features in their analysis order.
- `sn_input_feature_names.csv`: original H5 feature identifiers and symbols, including the duplicate-name mapping.
- `sn_factor_levels.csv`: categorical levels and their order.
- `GSE287257_donors.csv`: twelve public-study donors and their group, age and sex.

Script 3 constructs the count object using these fixed annotations and embedding. The preprocessing and annotation methods are described in the manuscript. Profiles require at least twenty nuclei per donor–cell-type combination, with at least three eligible donors per group for a cell type.

### Spatial transcriptomics

Download the eighteen Seurat objects listed in `config/GSE288365_samples.csv`. The `download_url` and `input_md5` columns identify each public file. This configuration also records donor, diagnosis, subtype, age, sex and eligibility for the comorbidity sensitivity analysis. Put only these eighteen `.rds.gz` objects in the spatial input directory.

Counts are combined within donor and supplied anatomical region. The four modeled regions are `Anterior_horns`, `Grey_matter`, `White_matter` and `Posterior_white_matter`. Region labels use `manual_3` when available and `manual` otherwise. The primary design is disease group. Sensitivity designs add sex, add age and sex, add standardized regional library size, or restrict the cohort by recorded comorbidity.

### Iron-related gene sets

`resources/iron_gene_sets.gmt` contains fourteen sets, 792 memberships and 545 unique symbols from Human MSigDB v2026.1.Hs. Overlap is retained. `source_data/Gene_set_label_key.csv` links full set names to S01–S14 and display labels. UCell uses maxRank=1,500, AUCell uses the top 5% of the expression universe, and singscore uses centered up-set scores.

## Run the transcriptomic analysis

```sh
Rscript --vanilla scripts/2_prepare_bulk.R "$ALS_INPUT/GSE287256_ALS_bulkRNAseq_results_bulk_counts_normalized.csv.gz" config/GSE287256_samples.csv "$ALS_OUT/bulk"
Rscript --vanilla scripts/3_prepare_single_nucleus.R "$ALS_INPUT/GSE287257_RAW" "$ALS_OUT/sn_input"
Rscript --vanilla scripts/4_aggregate_single_nucleus.R "$ALS_OUT/sn_input/GSE287257_annotated.rds" config/GSE287257_donors.csv "$ALS_OUT/sn_pseudobulk" "$ALS_OUT/sn_input/bpcells_raw"
Rscript --vanilla scripts/5_aggregate_spatial.R "$ALS_INPUT/GSE288365_objects" config/GSE288365_samples.csv "$ALS_OUT/spatial_counts"

Rscript --vanilla scripts/6_differential_expression.R "$ALS_OUT/bulk/bulk_expression.rds" "$ALS_OUT/sn_pseudobulk/sn_pseudobulk.rds" "$ALS_OUT/expression"
Rscript --vanilla scripts/7_single_nucleus_scores.R "$ALS_OUT/sn_input/GSE287257_annotated.rds" "$ALS_OUT/sn_pseudobulk/sn_pseudobulk.rds" "$ALS_OUT/nucleus_scores" "$ALS_OUT/sn_input/bpcells_raw"
Rscript --vanilla scripts/8_gene_set_models.R "$ALS_OUT/bulk/bulk_expression.rds" "$ALS_OUT/spatial_counts/spatial_pseudobulk.rds" config/GSE288365_samples.csv "$ALS_OUT/nucleus_scores/results/sn_donor_scores.csv" "$ALS_OUT/scores"
Rscript --vanilla scripts/9_spatial_expression.R "$ALS_OUT/spatial_counts/spatial_pseudobulk.rds" config/GSE288365_samples.csv "$ALS_OUT/expression/bulk_candidates.csv" "$ALS_OUT/spatial"
Rscript --vanilla scripts/10_spatial_lodo.R "$ALS_OUT/spatial_counts/spatial_pseudobulk.rds" config/GSE288365_samples.csv "$ALS_OUT/expression/bulk_candidates.csv" "$ALS_OUT/spatial/spatial_candidates_all_models.csv" "$ALS_OUT/spatial/spatial_supported_pairs.csv" "$ALS_OUT/lodo"
Rscript --vanilla scripts/11_spatial_qc_counts.R "$ALS_INPUT/GSE288365_objects" config/GSE288365_samples.csv "$ALS_OUT/spatial_counts/spatial_pseudobulk.rds" "$ALS_OUT/spatial_qc_counts"
Rscript --vanilla scripts/12_spatial_qc_models.R "$ALS_OUT/spatial_qc_counts" config/GSE288365_samples.csv "$ALS_OUT/expression/bulk_candidates.csv" "$ALS_OUT/spatial_qc_models"
Rscript --vanilla scripts/13_conditional_enrichment.R "$ALS_OUT/expression/bulk_DE_all.csv" "$ALS_OUT/spatial/spatial_candidates_all_models.csv" "$ALS_OUT/enrichment"
```

Script 7 scores individual nuclei before donor aggregation. Script 8 consumes those donor scores. LODO holds the candidate family and regional expression filters fixed while refitting each donor deletion. Enrichment uses one-sided Fisher tests, conditional odds ratios, two-sided exact 95% intervals and BH correction over the twenty-eight stage–set tests.

## Clinical associations

Use the authorized `ALS_Clinical_Data.xlsx`, sheet `Data1`, with unique nonmissing `Number` values. Patient records are subject to institutional approval and the manuscript's data-access conditions.

| Field | Column name or prefix | Coding |
|---|---|---|
| Patient key | `Number` | Unique patient identifier |
| Age | `Blood_collection_age` | Years |
| Sex | `Sex ` | F or M |
| Onset | `Onset_site` | Spinal, Bulbar or Mixture |
| Baseline score | `ALSFRS-R_baseline` | Numeric points |
| Follow-up score | `ALSFRS-R_follow_up` | Numeric points |
| Interval | `follow-up_interval` | Positive days |
| Ferritin | `Serum_ferritin (` | ng/mL |
| Serum iron | `serum_iron` | umol/L |
| TSAT | `transferrin_saturation` | Percent |
| Clinical remarks | `anemia`, `Inflammation`, `iron_supplements` | Blank/No: no recorded abnormality, Yes: recorded abnormality |
| Recorded slope | Header containing `ALSFRS-R/month` | Used to check the recomputed outcome's sign and value |

Monthly decline is `(baseline - follow-up) / days * 30.4`. Numeric follow-up scores and positive intervals define the evaluable outcome. Non-numeric missing tokens are missing values. The main analysis uses available marker/outcome pairs. Bootstrap intervals use 5,000 patient-paired resamples with seed 20260916. The strict sensitivity subset requires no recorded anemia, inflammation or iron supplementation. HC3 models adjust for age, sex, baseline ALSFRS-R and onset site.

```sh
Rscript --vanilla scripts/14_clinical_associations.R /authorized/path/ALS_Clinical_Data.xlsx "$ALS_OUT/clinical"
```

Outputs include `correlations.csv`, `adjusted_HC3.csv`, model diagnostics and the bootstrap specification. Patient-level rows are not exported.

## Figures

To render the manuscript's analysis figures directly from the included aggregate data:

```sh
Rscript --vanilla scripts/16_render_figures.R "$ALS_OUT/figures"
```

To render newly calculated results:

```sh
Rscript --vanilla scripts/15_prepare_plot_data.R "$ALS_OUT/expression" "$ALS_OUT/scores/score_tests.csv" "$ALS_OUT/spatial" "$ALS_OUT/enrichment" "$ALS_OUT/clinical" "$ALS_OUT/lodo/lodo" "$ALS_OUT/plot_data"
Rscript --vanilla scripts/16_render_figures.R "$ALS_OUT/figures_from_models" "$ALS_OUT/plot_data"
```

Script 15 uses explicit expression, score, spatial, enrichment, clinical and LODO outputs. UMAP, marker summaries, cell fractions, anatomical coordinates and display labels are fixed inputs. Each output figure is exported as PDF, SVG, PNG and TIFF.

The supplied `source_data/spatial_LODO/FigureS3_source_data.csv` is a compact gene–region–model summary for redrawing Figure S3. Given a full model-output directory, script 17 calculates that summary from `all_candidate_results.csv`. Its displayed ranges are minimum–maximum LODO effects, not confidence intervals.

```sh
Rscript --vanilla scripts/17_render_spatial_lodo.R "$ALS_OUT/lodo/lodo" "$ALS_OUT/figureS3"
```

## Resource attribution and access

GEO GSE287256, GSE287257 and GSE288365 remain available from their original public records. Cite the original studies and GEO accessions. Human MSigDB and its source collections retain their own [terms and citation requirements](https://www.gsea-msigdb.org/gsea/msigdb). The MIT license applies to original project code and does not replace third-party terms. Installed R libraries are obtained separately under their respective licenses. Clinical patient data follow the institutional access conditions stated in the manuscript.
