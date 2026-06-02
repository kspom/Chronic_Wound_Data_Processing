# R/DESeq2 replication of the differential expression analysis

Independent replication of the manuscript's DE analysis (originally implemented in MATLAB) using DESeq2 in R. Supports the methodological discussion in Supplementary Note S1 of the revised submission.

## Scripts (run in order)

`01_reconcile_samples.R` — sanity check: list Salmon `quant.sf` files and match against the sample metadata. Optional.

`02_import_originalT2G.R` — import Salmon quantifications, aggregate to gene-level counts, build the sample table, cache intermediates.

`03_DE_three_schemes.R` — run DESeq2 under three schemes: combined model with batch as covariate, per-batch + intersection, and pseudo-bulk per patient. Compare against the eight originally reported genes.

`04_qc_depth.R` — per-sample sequencing depth and gene-detection rates, stratified by batch, outcome, and treatment.

`05_treatment_covariate.R` — DE with treatment included as an additional covariate.

`06_sample_quality_figure.R` — 8-panel sample-quality figure (Figure S17): RIN, depth, and detected-gene relationships for both batches.

## Setup — required input files

Before running, place the following input files **in the same directory as the four R scripts**:

```
<your folder>/
├── 01_reconcile_samples.R
├── 02_import_originalT2G.R
├── 03_DE_three_schemes.R
├── 04_qc_depth.R
├── 05_treatment_covariate.R
├── 06_sample_quality_figure.R
├── README_DE_R.md
├── S4_Table_Chronic_Wound_Metadata.xlsx  ← copy from the data deposit
├── converted_transcripts_to_genes.csv    ← copy from the data deposit
├── S1_Table_Sample_Quality.xlsx          ← copy from the data deposit (used by 06_sample_quality_figure.R)
├── Cohort1/                              ← copy from the data deposit
│   ├── D11_CKDL220026941_quant.sf
│   ├── ...
│   └── (127 Salmon quant files)
└── Cohort2/                              ← copy from the data deposit
    ├── 659536_quant.sf
    ├── ...
    └── (309 Salmon quant files)
```

## Environment

R 4.4.3 with the following packages: `DESeq2` 1.46.0, `tximport` 1.34.0, `apeglm` 1.28.0, `dplyr`, `readxl`, `writexl`, `ggplot2`, `cowplot`. Installation via `conda`/`miniforge3` is recommended; all packages are available from the `conda-forge` and `bioconda` channels.
