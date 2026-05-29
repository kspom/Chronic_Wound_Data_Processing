# Build gene-level count matrix from Salmon quantifications.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tximport)
})

OUT <- "cache"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Tx-to-gene mapping
t2g <- read.csv("converted_transcripts_to_genes.csv", stringsAsFactors = FALSE)
cat("Mapping rows:", nrow(t2g),
    "  unique tx:", length(unique(t2g$ensembl_transcript_id)),
    "  unique genes:", length(unique(t2g$ensembl_gene_id)), "\n")
saveRDS(t2g, file.path(OUT, "tx2gene_originalT2G.rds"))
tx2gene <- t2g %>% select(ensembl_transcript_id, ensembl_gene_id)

# Sample manifest
meta <- read_excel("Chronic_Wound_Metadata.xlsx")
meta$healed <- gsub("Non-Healer", "Non-healer", meta$healed)

list_quants <- function(dir, id_strip) {
  f <- list.files(dir, pattern = "_quant\\.sf$", full.names = TRUE)
  f <- f[!grepl("/\\._", f)]
  setNames(f, sub(id_strip, "", basename(f)))
}
files_b1 <- list_quants("Cohort1", "_CKDL.*")
files_b2 <- list_quants("Cohort2", "_quant\\.sf$")
files_b1 <- files_b1[intersect(names(files_b1), meta$sample)]
files_b2 <- files_b2[intersect(names(files_b2), meta$sample)]
cat("Batch1 samples:", length(files_b1), "  Batch2:", length(files_b2), "\n")

# tximport per batch
import_batch <- function(files, label) {
  cache <- file.path(OUT, paste0("txi_originalT2G_", label, ".rds"))
  if (file.exists(cache)) { cat("Cached:", label, "\n"); return(readRDS(cache)) }
  cat("Importing", label, "...\n")
  txi <- tximport(files, type = "salmon",
                  tx2gene = tx2gene,
                  countsFromAbundance = "no",
                  ignoreAfterBar  = TRUE,
                  ignoreTxVersion = TRUE)
  saveRDS(txi, cache); txi
}
t1 <- import_batch(files_b1, "batch1"); gc()
t2 <- import_batch(files_b2, "batch2"); gc()

# Merge batches
stopifnot(all(rownames(t1$counts) == rownames(t2$counts)))
txi <- list(
  counts    = cbind(t1$counts,    t2$counts),
  abundance = cbind(t1$abundance, t2$abundance),
  length    = cbind(t1$length,    t2$length),
  countsFromAbundance = "no"
)
cat("Merged matrix:", nrow(txi$counts), "genes x", ncol(txi$counts), "samples\n")
saveRDS(txi, file.path(OUT, "txi_originalT2G_merged.rds"))

# Sample table
samp <- meta[match(colnames(txi$counts), meta$sample), ] %>%
  mutate(
    outcome   = factor(healed,    levels = c("Healer", "Non-healer")),
    batch     = factor(Batch,     levels = c("Batch1", "Batch2")),
    treatment = factor(Treatment, levels = c("Placebo", "Active", "StandCare")),
    week      = suppressWarnings(as.integer(sub("^w", "", week))),
    patient   = factor(patient)
  )
stopifnot(all(samp$sample == colnames(txi$counts)))
saveRDS(samp, file.path(OUT, "sample_table.rds"))

# Gene symbol map
gs <- t2g %>% distinct(ensembl_gene_id, external_gene_name)
gs <- data.frame(GENEID = gs$ensembl_gene_id,
                 SYMBOL = gs$external_gene_name,
                 stringsAsFactors = FALSE)
gs <- gs[match(rownames(txi$counts), gs$GENEID), ]
saveRDS(gs, file.path(OUT, "gene_symbols_originalT2G.rds"))

cat("\nDone. Cache written to:", OUT, "\n")
