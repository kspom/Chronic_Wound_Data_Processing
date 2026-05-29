# Figure S17: 8-panel RIN / depth / detected-genes scatterplots for both batches.
# Inputs:
#   cache/txi_originalT2G_*.rds (gene counts, built by 02_import_originalT2G.R)
#   Sample_Quality_Data.xlsx     (per-sample RIN + mapped reads, deposited as supplementary)
# Output:
#   results/Figure_S17.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(cowplot)
})

OUT <- "cache"
RES <- "results"
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

# Paper filter (zero in <= 150 of 436 samples) applied to merged counts
txi_m <- readRDS(file.path(OUT, "txi_originalT2G_merged.rds"))
counts_m <- round(txi_m$counts); storage.mode(counts_m) <- "integer"
keep_paper_merged <- rowSums(counts_m == 0) <= 150
genes_in_filter <- rownames(counts_m)[keep_paper_merged]
cat("Paper-filter gene set:", length(genes_in_filter), "genes\n")

# Per-batch detected-gene counts (raw and filter-passing)
gene_counts <- function(rds_path) {
  txi <- readRDS(rds_path)
  c_int <- round(txi$counts); storage.mode(c_int) <- "integer"
  c_filt <- c_int[rownames(c_int) %in% genes_in_filter, ]
  data.frame(sample = colnames(c_int),
             n_genes_raw  = colSums(c_int > 0),
             n_genes_filt = colSums(c_filt > 0))
}
gc1 <- gene_counts(file.path(OUT, "txi_originalT2G_batch1.rds"))
gc2 <- gene_counts(file.path(OUT, "txi_originalT2G_batch2.rds"))

# Sample quality data (per-sample RIN + mapped_reads_M)
sq <- read_excel("Sample_Quality_Data.xlsx")
sq <- sq %>% left_join(bind_rows(gc1, gc2), by = "sample")

b1 <- sq %>% filter(batch == "Batch1",
                    !is.na(RIN), !is.na(mapped_reads_M), !is.na(n_genes_raw))
b2 <- sq %>% filter(batch == "Batch2",
                    !is.na(RIN), !is.na(mapped_reads_M), !is.na(n_genes_raw))
cat(sprintf("N Batch 1 samples: %d   N Batch 2 samples: %d\n", nrow(b1), nrow(b2)))

cor_fmt <- function(x, y) {
  p <- cor.test(x, y, method = "pearson")
  s <- suppressWarnings(cor.test(x, y, method = "spearman"))
  sprintf("Pearson r=%+.2f (p=%.2g)\nSpearman rho=%+.2f (p=%.2g)\nn=%d",
          p$estimate, p$p.value, s$estimate, s$p.value, length(x))
}

# annot_pos: "tr" top-right, "tl" top-left, "br" bottom-right
make_panel <- function(df, xvar, yvar, title, xlab, ylab, color, annot_pos = "tr") {
  pos <- list(
    tr = c(x = Inf,  y = Inf,  hjust = 1.05,  vjust = 1.1),
    tl = c(x = -Inf, y = Inf,  hjust = -0.05, vjust = 1.1),
    br = c(x = Inf,  y = -Inf, hjust = 1.05,  vjust = -0.05)
  )[[annot_pos]]
  ggplot(df, aes(.data[[xvar]], .data[[yvar]])) +
    geom_point(alpha = 0.5, size = 1.2, color = color) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.4) +
    annotate("text", x = pos["x"], y = pos["y"],
             hjust = pos["hjust"], vjust = pos["vjust"], size = 3.1,
             label = cor_fmt(df[[xvar]], df[[yvar]])) +
    labs(title = title, x = xlab, y = ylab) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 11, hjust = 0),
          aspect.ratio = 1 / 1.7,
          plot.margin = margin(2, 2, 2, 2))
}

b1$n_genes_raw_k  <- b1$n_genes_raw / 1000
b1$n_genes_filt_k <- b1$n_genes_filt / 1000
b2$n_genes_raw_k  <- b2$n_genes_raw / 1000
b2$n_genes_filt_k <- b2$n_genes_filt / 1000

C1 <- "steelblue"; C2 <- "tomato"

p1A  <- make_panel(b1, "RIN", "n_genes_raw_k",
                   "A. Batch 1: RIN vs detected genes (raw)",
                   "RIN", "Detected genes (k)", C1, "tr")
p2A  <- make_panel(b2, "RIN", "n_genes_raw_k",
                   "B. Batch 2: RIN vs detected genes (raw)",
                   "RIN", "Detected genes (k)", C2, "br")
p1Ap <- make_panel(b1, "RIN", "n_genes_filt_k",
                   "C. Batch 1: RIN vs genes (filtered 18,577 set)",
                   "RIN", "Detected genes (k)", C1, "br")
p2Ap <- make_panel(b2, "RIN", "n_genes_filt_k",
                   "D. Batch 2: RIN vs genes (filtered 18,577 set)",
                   "RIN", "Detected genes (k)", C2, "br")
p1B  <- make_panel(b1, "RIN", "mapped_reads_M",
                   "E. Batch 1: RIN vs sequencing depth",
                   "RIN", "Mapped reads (M)", C1, "tr")
p2B  <- make_panel(b2, "RIN", "mapped_reads_M",
                   "F. Batch 2: RIN vs sequencing depth",
                   "RIN", "Mapped reads (M)", C2, "tl")
p1C  <- make_panel(b1, "mapped_reads_M", "n_genes_raw_k",
                   "G. Batch 1: depth vs detected genes",
                   "Mapped reads (M)", "Detected genes (k)", C1, "tr")
p2C  <- make_panel(b2, "mapped_reads_M", "n_genes_raw_k",
                   "H. Batch 2: depth vs detected genes",
                   "Mapped reads (M)", "Detected genes (k)", C2, "br")

merged <- plot_grid(
  p1A,  p2A,
  p1Ap, p2Ap,
  p1B,  p2B,
  p1C,  p2C,
  ncol = 2, align = "hv"
)
ggsave(file.path(RES, "Figure_S17.pdf"), merged, width = 9, height = 11)
cat("Written:", file.path(RES, "Figure_S17.pdf"), "\n")
