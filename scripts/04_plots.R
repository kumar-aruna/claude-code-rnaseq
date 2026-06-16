# ============================================================
# Step 4: Visualization
# Goal:
#   1. Volcano plot for each of the 5 comparisons
#      (x = log2 fold-change, y = -log10 adjusted p-value)
#   2. Heatmap of the top hypoxia-responsive genes across samples
# Output: PNG figures in results/
# ============================================================

suppressMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)   # auto-spaces overlapping gene labels
})

# --- A reusable volcano-plot function --------------------------------
# A volcano plot shows EVERY gene as a dot:
#   - far left/right  = big fold-change (down / up)
#   - high up         = very statistically significant
# So the interesting genes are in the top corners.
volcano <- function(csv, title, out_png) {
  d <- read.csv(csv, row.names = 1)
  d <- d[!is.na(d$padj), ]                       # drop genes with no padj

  # Some very strong genes have padj that underflows to exactly 0, which
  # would make -log10(padj) = Inf and silently drop them from the plot.
  # Floor those at the smallest non-zero padj so they stay visible.
  if (any(d$padj == 0)) d$padj[d$padj == 0] <- min(d$padj[d$padj > 0])

  # classify each gene for coloring
  d$status <- "Not sig"
  d$status[d$padj < 0.05 & d$log2FoldChange >  1] <- "Up"
  d$status[d$padj < 0.05 & d$log2FoldChange < -1] <- "Down"

  # label the 10 most significant genes
  top <- head(d[order(d$padj), ], 10)
  top$gene <- rownames(top)

  p <- ggplot(d, aes(log2FoldChange, -log10(padj), color = status)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("Up" = "firebrick",
                                  "Down" = "steelblue",
                                  "Not sig" = "grey75")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_text_repel(data = top, aes(label = gene), size = 3,
                    color = "black", max.overlaps = Inf,
                    box.padding = 0.5, min.segment.length = 0) +
    labs(title = title, x = "log2 fold-change", y = "-log10(adjusted p)") +
    theme_bw() +
    theme(legend.title = element_blank())

  ggsave(out_png, p, width = 7, height = 5.5, dpi = 150)
  cat("Saved", out_png,
      sprintf("(Up: %d, Down: %d)\n",
              sum(d$status == "Up"), sum(d$status == "Down")))
}

# --- Build a volcano for each comparison -----------------------------
if (!file.exists("results/DE_1_hypoxia_vs_normoxia_in_control.csv"))
  stop("Missing DE result CSVs - run scripts/03_deseq.R first.")

volcano("results/DE_1_hypoxia_vs_normoxia_in_control.csv",
        "Hypoxia vs Normoxia (Control)", "results/volcano_1_hypoxia.png")
volcano("results/DE_2_ITPR3KO_vs_control.csv",
        "ITPR3 KO vs Control (Normoxia)", "results/volcano_2_ITPR3KO.png")
volcano("results/DE_3_RELBKO_vs_control.csv",
        "RELB KO vs Control (Normoxia)", "results/volcano_3_RELBKO.png")
volcano("results/DE_4_ITPR3_hypoxia_interaction.csv",
        "ITPR3 x Hypoxia interaction", "results/volcano_4_ITPR3_interaction.png")
volcano("results/DE_5_RELB_hypoxia_interaction.csv",
        "RELB x Hypoxia interaction", "results/volcano_5_RELB_interaction.png")

# --- Heatmap of top hypoxia genes ------------------------------------
# Take the 30 most significant genes from comparison #1 and show their
# (scaled) expression across all 12 samples. Replicates should look alike
# and hypoxia vs normoxia should form two clear blocks.
dds <- readRDS("results/dds_step1.rds")
vsd <- vst(dds, blind = TRUE)

# Read the DE table once, order by padj, take the top 30 gene names.
de1 <- read.csv("results/DE_1_hypoxia_vs_normoxia_in_control.csv", row.names = 1)
de1 <- de1[order(de1$padj), ]
top_genes <- head(rownames(de1), 30)
top_genes <- intersect(top_genes, rownames(vsd))   # keep only genes present in vsd

mat <- assay(vsd)[top_genes, ]
mat <- t(scale(t(mat)))                          # z-score each gene (row)

ann <- as.data.frame(colData(vsd)[, c("genotype", "oxygen")])

png("results/heatmap_top_hypoxia_genes.png", width = 900, height = 950, res = 120)
pheatmap(mat,
         annotation_col = ann,
         show_colnames = FALSE,
         main = "Top 30 hypoxia-responsive genes (z-scored)")
dev.off()
cat("Saved results/heatmap_top_hypoxia_genes.png\n")

cat("\nStep 4 complete.\n")
