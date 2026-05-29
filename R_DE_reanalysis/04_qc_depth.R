# Sequencing depth (Salmon-mapped reads) per sample, summarized by batch / outcome / treatment.

suppressPackageStartupMessages({
  library(dplyr)
  library(writexl)
})

OUT <- "cache"
RES <- "results"
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

txi  <- readRDS(file.path(OUT, "txi_originalT2G_merged.rds"))
samp <- readRDS(file.path(OUT, "sample_table.rds"))

samp$mapped_reads   <- colSums(txi$counts)
samp$detected_genes <- colSums(txi$counts >= 1)

# Per-sample table
per_sample <- samp %>%
  select(sample, patient, batch, outcome, treatment,
         mapped_reads, detected_genes) %>%
  mutate(mapped_reads_M = round(mapped_reads / 1e6, 2))

# Group summaries
summarize_depth <- function(df, group_vars) {
  df %>% group_by(across(all_of(group_vars))) %>%
    summarise(
      n_samples             = n(),
      median_reads_M        = round(median(mapped_reads)/1e6, 2),
      mean_reads_M          = round(mean(mapped_reads)/1e6, 2),
      sd_reads_M            = round(sd(mapped_reads)/1e6, 2),
      min_reads_M           = round(min(mapped_reads)/1e6, 2),
      max_reads_M           = round(max(mapped_reads)/1e6, 2),
      median_genes_detected = median(detected_genes),
      .groups = "drop"
    )
}
by_batch                   <- summarize_depth(samp, "batch")
by_batch_outcome           <- summarize_depth(samp, c("batch","outcome"))
by_batch_treatment         <- summarize_depth(samp, c("batch","treatment"))
by_batch_outcome_treatment <- summarize_depth(samp, c("batch","outcome","treatment"))

cat("=== By batch ===\n");                          print(by_batch)
cat("\n=== By batch x outcome ===\n");              print(by_batch_outcome)
cat("\n=== By batch x treatment ===\n");            print(by_batch_treatment)
cat("\n=== By batch x outcome x treatment ===\n");  print(by_batch_outcome_treatment)

# Wilcoxon test: depth difference healer vs non-healer
cat("\n=== Wilcoxon: mapped reads, healer vs non-healer, per batch ===\n")
for (b in levels(samp$batch)) {
  sub <- samp[samp$batch == b, ]
  if (dplyr::n_distinct(sub$outcome) < 2) next
  w <- wilcox.test(mapped_reads ~ outcome, data = sub)
  cat(sprintf("  %s:  W=%g,  p=%.4f\n", b, w$statistic, w$p.value))
}

# Save
write_xlsx(
  list(
    per_sample                 = per_sample,
    summary_by_batch           = by_batch,
    summary_by_batch_outcome   = by_batch_outcome,
    summary_by_batch_treatment = by_batch_treatment,
    summary_by_b_o_t           = by_batch_outcome_treatment
  ),
  file.path(RES, "QC_sequencing_depth.xlsx")
)
cat("\nWritten:", file.path(RES, "QC_sequencing_depth.xlsx"), "\n")
