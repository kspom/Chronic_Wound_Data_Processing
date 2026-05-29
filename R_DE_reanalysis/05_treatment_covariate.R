# Differential expression with treatment included as a model covariate.
# Combined per-sample (~ batch + treatment + outcome) and pseudo-bulk per patient.

suppressPackageStartupMessages({
  library(dplyr)
  library(DESeq2)
  library(writexl)
})

OUT <- "cache"
RES <- "results"
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

txi  <- readRDS(file.path(OUT, "txi_originalT2G_merged.rds"))
samp <- readRDS(file.path(OUT, "sample_table.rds"))
gs   <- readRDS(file.path(OUT, "gene_symbols_originalT2G.rds"))

# Gene filter (same as paper)
counts_int <- round(txi$counts); storage.mode(counts_int) <- "integer"
keep_paper <- rowSums(counts_int == 0) <= 150
counts_filt <- counts_int[keep_paper, ]
abund_filt  <- txi$abundance[keep_paper, ]
len_filt    <- txi$length[keep_paper, ]
gs_filt     <- gs[keep_paper, ]
rownames(counts_filt) <- gs_filt$GENEID
rownames(abund_filt)  <- gs_filt$GENEID
rownames(len_filt)    <- gs_filt$GENEID

report_de <- function(res, label, cuts) {
  cat(sprintf("\n  %s:\n", label))
  for (k in seq_len(nrow(cuts))) {
    p <- cuts$padj[k]; lf <- cuts$lfc[k]
    ok <- !is.na(res$padj) & res$padj < p & abs(res$log2FoldChange) >= lf
    cat(sprintf("    padj<%.2f & |lfc|>=%g:  %5d  (%5d up, %5d down)\n",
                p, lf, sum(ok),
                sum(ok & res$log2FoldChange > 0),
                sum(ok & res$log2FoldChange < 0)))
  }
}
cuts <- data.frame(padj = c(0.01, 0.01, 0.05, 0.05),
                   lfc  = c(0,    1,    0,    1))

# Combined per-sample
cat("\n=== Combined per-sample, ~ batch + treatment + outcome ===\n")
ddsAt <- DESeqDataSetFromTximport(
  list(counts = counts_filt, abundance = abund_filt, length = len_filt,
       countsFromAbundance = "no"),
  colData = samp, design = ~ batch + treatment + outcome
)
ddsAt <- DESeq(ddsAt)
resAt <- results(ddsAt, name = "outcome_Non.healer_vs_Healer", alpha = 0.05)
report_de(resAt, "combined +treatment", cuts)

# Pseudo-bulk per patient
cat("\n=== Pseudo-bulk per patient, ~ batch + treatment + outcome ===\n")
samp$pb_id <- paste(samp$patient, samp$batch, sep = "::")
pb_groups <- split(seq_len(nrow(samp)), samp$pb_id)
pb_counts <- sapply(pb_groups, function(idx) rowSums(counts_filt[, idx, drop = FALSE]))
storage.mode(pb_counts) <- "integer"

get_modal <- function(x) names(sort(table(x), decreasing = TRUE))[1]
pb_samp <- samp %>% group_by(pb_id) %>%
  summarise(patient   = dplyr::first(patient),
            batch     = dplyr::first(batch),
            outcome   = dplyr::first(outcome),
            treatment = factor(get_modal(as.character(treatment)),
                               levels = levels(samp$treatment)),
            .groups   = "drop")
pb_samp <- pb_samp[match(colnames(pb_counts), pb_samp$pb_id), ]

ddsCt <- DESeqDataSetFromMatrix(pb_counts, colData = pb_samp,
                                design = ~ batch + treatment + outcome)
ddsCt <- DESeq(ddsCt)
resCt <- results(ddsCt, name = "outcome_Non.healer_vs_Healer", alpha = 0.05)
report_de(resCt, "pseudo-bulk +treatment", cuts)

# Treatment effect (sanity check)
cat("\n=== Treatment effect (Active vs Placebo) ===\n")
for (nm in resultsNames(ddsAt)) {
  if (grepl("treatment", nm)) {
    r <- results(ddsAt, name = nm, alpha = 0.05)
    ok <- !is.na(r$padj) & r$padj < 0.05
    cat(sprintf("  %s:  %d genes at padj<0.05\n", nm, sum(ok)))
  }
}

# Save
attach_symbol <- function(res) {
  d <- as.data.frame(res)
  d$gene_id <- rownames(d)
  d$SYMBOL  <- gs$SYMBOL[match(d$gene_id, gs$GENEID)]
  d %>% select(gene_id, SYMBOL, baseMean, log2FoldChange, padj, stat)
}
write_xlsx(
  list(
    combined_treatment   = attach_symbol(resAt)[order(attach_symbol(resAt)$padj), ],
    pseudobulk_treatment = attach_symbol(resCt)[order(attach_symbol(resCt)$padj), ]
  ),
  file.path(RES, "DE_treatment_covariate.xlsx")
)

cat("\nDone.\n")
