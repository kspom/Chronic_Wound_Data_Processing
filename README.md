# Chronic_Wound_Data_Processing
Analysis of longitudinal transcriptomic data from human chronic wounds (diabetic foot ulcers) collected on a weekly basis.

## Data reading and analysis (run in this order)

| Script | Purpose |
|---|---|
| `Read_Data` | reads Salmon quantifications and metadata; output `ChronicWoundData.mat` |
| `BoxPlotSomeGenes` | a tool for plotting desired genes |
| `BatchEffectCorrection` | log-transformation and batch-wise standardisation; output `ChronicWoundDataBatchCorrected.mat` |
| `Clusters_analysis` | deriving clusters, plotting scatterplots and cluster time-series |
| `Timeseries_plot` | time series of gene expression, non-normalized, one cluster, each patient |
| `validation_pcr_scatters` | comparing RNAseq and qPCR, scatterplots |
| `qpcr_validation_barplots` | comparing RNAseq and qPCR, barplots |
| `PrinCompAnalysis` | PCA for non-normalized gene expression data |
| `markov_chain_model` | Markov Chain Model for healers and non-healers |
| `DEanalysis` | differentially expressed genes between healers and non-healers: for all, treated and non-treated samples |

## Sensitivity analyses (Fig S24, S15–S16 Tables)

| Script | Purpose |
|---|---|
| `PatientResampling` | patient-level bootstrap, permutation test and leave-one-patient-out for transition probabilities and group differences; input `Patient_transitions.xlsx` (output of `markov_chain_model`), output `Resampling_results.xlsx` |
| `RemoveLowRIN` | removes samples with RIN below a cut-off (`RINcut`); input `ChronicWoundData.mat` and `S1_Table_Sample_Quality.xlsx`, output `ChronicWoundData_lowRIN.mat` |
| `NormalizeCounts` | median-of-ratios library-size normalization of counts; input `ChronicWoundData.mat`, output `ChronicWoundData_normalized.mat` |

To run the cluster and Markov-chain analysis on either dataset, switch the `load` line at the top of `BatchEffectCorrection` to the corresponding file (commented lines are provided), then run `BatchEffectCorrection` and `markov_chain_model`.

## R re-analysis of differential expression

The folder `R_DE_reanalysis` contains the DESeq2-based replication of the differential expression analysis (S1 Text); see `R_DE_reanalysis/README_DE_R.md`.
