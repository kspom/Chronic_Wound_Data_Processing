# Sanity check: match Salmon quant.sf files to metadata sample IDs.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

m <- read_excel("Chronic_Wound_Metadata.xlsx")
m$healed <- gsub("Non-Healer", "Non-healer", m$healed)

# Batch 1
b1files <- list.files("Cohort1", pattern = "_quant\\.sf$")
b1files <- b1files[!grepl("^\\._", b1files)]
b1ids   <- sub("_CKDL.*", "", b1files)

cat("Batch1 folder samples:", length(b1ids), "\n")
b1meta <- m %>% filter(Batch == "Batch1") %>% pull(sample)
cat("Batch1 metadata samples:", length(b1meta), "\n")
cat("In folder but not metadata:\n"); print(setdiff(b1ids, b1meta))
cat("In metadata but not folder:\n"); print(setdiff(b1meta, b1ids))

# Batch 2
b2files <- list.files("Cohort2", pattern = "_quant\\.sf$")
b2files <- b2files[!grepl("^\\._", b2files)]
b2ids   <- sub("_quant\\.sf$", "", b2files)

cat("\nBatch2 folder samples:", length(b2ids), "\n")
b2meta <- m %>% filter(Batch == "Batch2") %>% pull(sample)
cat("Batch2 metadata samples:", length(b2meta), "\n")
cat("In folder but not metadata:\n"); print(setdiff(b2ids, b2meta))
cat("In metadata but not folder:\n"); print(setdiff(b2meta, b2ids))
