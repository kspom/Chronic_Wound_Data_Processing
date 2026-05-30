# Differential expression: three schemes (combined, per-batch, pseudo-bulk).

suppressPackageStartupMessages({
  library(dplyr)
  library(DESeq2)
  library(writexl)
})

OUT  <- "cache"
RES  <- "results"
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

txi  <- readRDS(file.path(OUT, "txi_originalT2G_merged.rds"))
samp <- readRDS(file.path(OUT, "sample_table.rds"))
gs   <- readRDS(file.path(OUT, "gene_symbols_originalT2G.rds"))

# Gene filter
counts_int <- round(txi$counts); storage.mode(counts_int) <- "integer"
zero_per_gene <- rowSums(counts_int == 0)
keep_paper <- zero_per_gene <= 150
cat("Genes kept:", sum(keep_paper), "  (paper reports 18,924)\n\n")

counts_filt <- counts_int[keep_paper, ]
abund_filt  <- txi$abundance[keep_paper, ]
len_filt    <- txi$length[keep_paper, ]
gs_filt     <- gs[keep_paper, ]
rownames(counts_filt) <- gs_filt$GENEID
rownames(abund_filt)  <- gs_filt$GENEID
rownames(len_filt)    <- gs_filt$GENEID

# Thresholds
report_de <- function(res, label, cuts) {
  cat(sprintf("\n  %s:\n", label))
  for (k in seq_len(nrow(cuts))) {
    p <- cuts$padj[k]; lf <- cuts$lfc[k]
    ok <- !is.na(res$padj) & res$padj < p & abs(res$log2FoldChange) >= lf
    up <- ok & res$log2FoldChange > 0
    dn <- ok & res$log2FoldChange < 0
    cat(sprintf("    padj<%.2f & |lfc|>=%g:  %5d  (%5d up, %5d down)\n",
                p, lf, sum(ok), sum(up), sum(dn)))
  }
}
cuts <- data.frame(padj = c(0.01, 0.01, 0.05, 0.05),
                   lfc  = c(0,    1,    0,    1))

# Scheme A: combined per-sample
cat("\n=== A. Combined ~ batch + outcome ===\n")
ddsA <- DESeqDataSetFromTximport(
  list(counts = counts_filt, abundance = abund_filt, length = len_filt,
       countsFromAbundance = "no"),
  colData = samp, design = ~ batch + outcome
)
ddsA <- DESeq(ddsA)
resA <- results(ddsA, name = "outcome_Non.healer_vs_Healer", alpha = 0.05)
report_de(resA, "combined", cuts)

# Scheme B: per-batch + intersection
cat("\n=== B. Per-batch + intersection ===\n")
run_one_batch <- function(b) {
  idx <- which(samp$batch == b)
  dds <- DESeqDataSetFromTximport(
    list(counts = counts_filt[, idx], abundance = abund_filt[, idx],
         length = len_filt[, idx], countsFromAbundance = "no"),
    colData = samp[idx, ], design = ~ outcome
  )
  dds <- DESeq(dds)
  results(dds, name = "outcome_Non.healer_vs_Healer", alpha = 0.05)
}
res1 <- run_one_batch("Batch1")
res2 <- run_one_batch("Batch2")
report_de(res1, "Batch1", cuts)
report_de(res2, "Batch2", cuts)

cat("\n  Intersection:\n")
for (k in seq_len(nrow(cuts))) {
  p <- cuts$padj[k]; lf <- cuts$lfc[k]
  up1 <- rownames(res1)[!is.na(res1$padj) & res1$padj < p & res1$log2FoldChange >=  lf]
  dn1 <- rownames(res1)[!is.na(res1$padj) & res1$padj < p & res1$log2FoldChange <= -lf]
  up2 <- rownames(res2)[!is.na(res2$padj) & res2$padj < p & res2$log2FoldChange >=  lf]
  dn2 <- rownames(res2)[!is.na(res2$padj) & res2$padj < p & res2$log2FoldChange <= -lf]
  ov_up <- length(intersect(up1, up2))
  ov_dn <- length(intersect(dn1, dn2))
  cat(sprintf("    padj<%.2f & |lfc|>=%g:  total %d  (up %d, down %d)\n",
              p, lf, ov_up + ov_dn, ov_up, ov_dn))
}

# Scheme C: pseudo-bulk per patient
cat("\n=== C. Pseudo-bulk per patient ===\n")
samp$pb_id <- paste(samp$patient, samp$batch, sep = "::")
pb_groups <- split(seq_len(nrow(samp)), samp$pb_id)
pb_counts <- sapply(pb_groups, function(idx) rowSums(counts_filt[, idx, drop = FALSE]))
storage.mode(pb_counts) <- "integer"
pb_samp <- samp %>% group_by(pb_id) %>%
  summarise(patient = dplyr::first(patient),
            batch   = dplyr::first(batch),
            outcome = dplyr::first(outcome),
            .groups = "drop")
pb_samp <- pb_samp[match(colnames(pb_counts), pb_samp$pb_id), ]

ddsC <- DESeqDataSetFromMatrix(pb_counts, colData = pb_samp,
                               design = ~ batch + outcome)
ddsC <- DESeq(ddsC)
resC <- results(ddsC, name = "outcome_Non.healer_vs_Healer", alpha = 0.05)
report_de(resC, "pseudo-bulk", cuts)

# Original 8 genes
orig_up <- c("CRTAC1")
orig_dn <- c("FLG","FLG2","IL36A","SPRR2F","LORICRIN","CALML5","SERPINB4")

attach_symbol <- function(res) {
  d <- as.data.frame(res)
  d$gene_id <- rownames(d)
  d$SYMBOL  <- gs$SYMBOL[match(d$gene_id, gs$GENEID)]
  d
}
dA  <- attach_symbol(resA)
dC  <- attach_symbol(resC)
d1  <- attach_symbol(res1)
d2  <- attach_symbol(res2)

per_batch_int <- function(res1, res2, padj_cut, lfc_cut) {
  up1 <- rownames(res1)[!is.na(res1$padj) & res1$padj < padj_cut & res1$log2FoldChange >=  lfc_cut]
  dn1 <- rownames(res1)[!is.na(res1$padj) & res1$padj < padj_cut & res1$log2FoldChange <= -lfc_cut]
  up2 <- rownames(res2)[!is.na(res2$padj) & res2$padj < padj_cut & res2$log2FoldChange >=  lfc_cut]
  dn2 <- rownames(res2)[!is.na(res2$padj) & res2$padj < padj_cut & res2$log2FoldChange <= -lfc_cut]
  list(up = intersect(up1, up2), dn = intersect(dn1, dn2))
}
pbi <- per_batch_int(res1, res2, 0.01, 0)
pbi_up <- gs$SYMBOL[match(pbi$up, gs$GENEID)]
pbi_dn <- gs$SYMBOL[match(pbi$dn, gs$GENEID)]

cat("\n=== Original 8 genes in each analysis ===\n")
show_gene <- function(d, sym) {
  r <- d[d$SYMBOL == sym, ]
  if (nrow(r) == 0) return(c(NA, NA))
  c(r$log2FoldChange[1], r$padj[1])
}
orig_genes <- c(orig_up, orig_dn)
tab <- data.frame(
  SYMBOL        = orig_genes,
  Combined_lfc  = sapply(orig_genes, function(g) round(show_gene(dA, g)[1], 2)),
  Combined_padj = sapply(orig_genes, function(g) signif(show_gene(dA, g)[2], 3)),
  PB_lfc        = sapply(orig_genes, function(g) round(show_gene(dC, g)[1], 2)),
  PB_padj       = sapply(orig_genes, function(g) signif(show_gene(dC, g)[2], 3)),
  In_per_batch_intersection = ifelse(orig_genes %in% c(pbi_up, pbi_dn), "yes", "no")
)
print(tab, row.names = FALSE)

cat("\n=== Per-batch + intersection gene list ===\n")
cat("  Up:  ", paste(pbi_up, collapse = ", "), "\n")
cat("  Down:", paste(pbi_dn, collapse = ", "), "\n")

# Save
saveRDS(list(resA = resA, res1 = res1, res2 = res2, resC = resC,
             gs = gs, samp = samp),
        file.path(OUT, "originalT2G_DE_final.rds"))

cols <- c("gene_id","SYMBOL","baseMean","log2FoldChange","padj","stat")

# Per-batch + intersection (DESeq2 replication of original MATLAB scheme)
intersection_tbl <- data.frame(
  gene_id   = c(pbi$up, pbi$dn),
  SYMBOL    = c(pbi_up, pbi_dn),
  direction = c(rep("up_in_non.healer", length(pbi$up)),
                rep("down_in_non.healer", length(pbi$dn))),
  Batch1_lfc  = c(res1$log2FoldChange[match(pbi$up, rownames(res1))],
                  res1$log2FoldChange[match(pbi$dn, rownames(res1))]),
  Batch1_padj = c(res1$padj[match(pbi$up, rownames(res1))],
                  res1$padj[match(pbi$dn, rownames(res1))]),
  Batch2_lfc  = c(res2$log2FoldChange[match(pbi$up, rownames(res2))],
                  res2$log2FoldChange[match(pbi$dn, rownames(res2))]),
  Batch2_padj = c(res2$padj[match(pbi$up, rownames(res2))],
                  res2$padj[match(pbi$dn, rownames(res2))])
)

write_xlsx(
  list(
    combined_per_sample    = dA[order(dA$padj), cols],
    per_batch_batch1       = d1[order(d1$padj), cols],
    per_batch_batch2       = d2[order(d2$padj), cols],
    per_batch_intersection = intersection_tbl,
    pseudobulk             = dC[order(dC$padj), cols],
    original_8_in_new      = tab
  ),
  file.path(RES, "DE_originalT2G_final.xlsx")
)

cat("\n=== Per-batch + intersection table (saved as sheet per_batch_intersection) ===\n")
print(intersection_tbl, row.names = FALSE)

cat("\nDone.\n")
