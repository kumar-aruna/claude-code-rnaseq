# ============================================================
# Step 2: Quality control (QC)
# Goal:
#   1. Transform the counts so they are suitable for plotting
#   2. PCA plot  -> do samples cluster by condition? (did the experiment work?)
#   3. Sample-to-sample distance heatmap -> spot any outlier/swapped samples
# Output: figures saved into results/
# ============================================================

suppressMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

# --- Load the object we built in Step 1 ---
if (!file.exists("results/dds_step1.rds"))
  stop("Missing results/dds_step1.rds - run scripts/01_setup.R first.")
dds <- readRDS("results/dds_step1.rds")

# --- 1. Variance-stabilizing transform (VST) -------------------------
# Raw counts are NOT good for distance/PCA plots: high-count genes would
# dominate just because their numbers are big. VST puts all genes on a
# comparable scale so the plot reflects real biological differences.
vsd <- vst(dds, blind = TRUE)   # blind = ignore the design, pure QC view

# --- 2. PCA plot -----------------------------------------------------
# PCA squashes 20,759 genes down to 2 summary axes (PC1, PC2). Samples
# that behave similarly land near each other. We expect samples to split
# by oxygen and/or genotype -- that means the biology is the main signal.
pca <- plotPCA(vsd, intgroup = c("genotype", "oxygen"), returnData = TRUE)
pct <- round(100 * attr(pca, "percentVar"))   # variance explained per axis

p_pca <- ggplot(pca, aes(PC1, PC2, color = oxygen, shape = genotype)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", pct[1], "% variance")) +
  ylab(paste0("PC2: ", pct[2], "% variance")) +
  ggtitle("PCA of GSE197576 samples") +
  theme_bw()

ggsave("results/qc_pca.png", p_pca, width = 7, height = 5, dpi = 150)
cat("Saved results/qc_pca.png\n")

# --- 3. Sample-to-sample distance heatmap ----------------------------
# Compute how different each sample is from every other sample.
# Replicates of the same condition should be the most similar (darkest).
sample_dists <- dist(t(assay(vsd)))
dist_mat <- as.matrix(sample_dists)

png("results/qc_sample_distances.png", width = 800, height = 700, res = 120)
pheatmap(dist_mat,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         main = "Sample-to-sample distances")
dev.off()
cat("Saved results/qc_sample_distances.png\n")

# --- Quick text summary of PC1 vs condition --------------------------
cat("\nPCA coordinates (for a quick sanity read):\n")
print(pca[, c("oxygen", "genotype", "PC1", "PC2")])

cat("\nStep 2 complete.\n")
