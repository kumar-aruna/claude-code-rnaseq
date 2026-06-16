# ============================================================
# Step 3: Differential expression (the core result)
# Goal: run the DESeq2 model and extract the 5 biologically
#       meaningful comparisons, saving each as a CSV table.
#
# Reminder on what a "comparison" is: we pick two groups of samples,
# and DESeq2 tests EVERY gene for whether it differs between them.
# Output = a ranked list of genes (we don't choose the genes).
# ============================================================

suppressMessages(library(DESeq2))

if (!file.exists("results/dds_step1.rds"))
  stop("Missing results/dds_step1.rds - run scripts/01_setup.R first.")
dds <- readRDS("results/dds_step1.rds")

# --- Run the model ---------------------------------------------------
# DESeq() estimates, for every gene, how its expression depends on
# genotype, oxygen, and their interaction (design = ~ genotype * oxygen).
dds <- DESeq(dds)

# The model's coefficients (the building blocks of our comparisons):
cat("Model coefficients available:\n")
print(resultsNames(dds))
cat("\n")

# --- Helper: extract one comparison, save ----------------------------
# NOTE on shrinkage: we report UNSHRUNKEN log2 fold-changes here. The
# 'normal' shrinkage method does not support interaction designs, and
# apeglm/ashr aren't installed on this machine. Shrinkage only adjusts
# the fold-change magnitudes of low-count genes for nicer ranking/plots;
# it does NOT change p-values or which genes are significant. This is a
# valid, commonly-reported analysis. (Optional future polish: install
# apeglm and re-run with lfcShrink.)
extract <- function(coef_name, out_csv, question) {
  res <- results(dds, name = coef_name)
  res <- res[order(res$padj), ]                      # best (smallest padj) first
  write.csv(as.data.frame(res), out_csv)
  n_sig <- sum(res$padj < 0.05, na.rm = TRUE)        # significant genes (FDR 5%)
  cat(sprintf("%-45s | sig genes (padj<0.05): %5d -> %s\n",
              question, n_sig, out_csv))
  invisible(res)
}

cat("Running the 5 comparisons:\n")
cat("------------------------------------------------------------\n")

# 1. Hypoxia vs Normoxia in CONTROL cells (the core hypoxia response).
#    In a "~ genotype * oxygen" model, the plain oxygen coefficient is
#    the hypoxia effect in the reference genotype (Control).
res_hyp <- extract("oxygen_Hypoxia_vs_Normoxia",
                   "results/DE_1_hypoxia_vs_normoxia_in_control.csv",
                   "1. Hypoxia vs Normoxia (Control)")

# 2. ITPR3 knockout vs Control (under normoxia baseline).
extract("genotype_ITPR3_KO_vs_Control",
        "results/DE_2_ITPR3KO_vs_control.csv",
        "2. ITPR3_KO vs Control (Normoxia)")

# 3. RELB knockout vs Control (under normoxia baseline).
extract("genotype_RELB_KO_vs_Control",
        "results/DE_3_RELBKO_vs_control.csv",
        "3. RELB_KO vs Control (Normoxia)")

# 4a. INTERACTION: does losing ITPR3 change the hypoxia response?
#     Genes here respond to hypoxia DIFFERENTLY in ITPR3_KO vs Control.
extract("genotypeITPR3_KO.oxygenHypoxia",
        "results/DE_4_ITPR3_hypoxia_interaction.csv",
        "4a. ITPR3 x hypoxia interaction")

# 4b. INTERACTION: does losing RELB change the hypoxia response?
extract("genotypeRELB_KO.oxygenHypoxia",
        "results/DE_5_RELB_hypoxia_interaction.csv",
        "4b. RELB x hypoxia interaction")

cat("------------------------------------------------------------\n\n")

# --- Biological sanity check ----------------------------------------
# If our pipeline is correct, canonical hypoxia genes must be strongly
# UP in hypoxia (positive log2FoldChange, tiny padj) in comparison #1.
cat("Sanity check - known hypoxia marker genes in comparison #1:\n")
markers <- c("VEGFA", "CA9", "SLC2A1", "PGK1", "LDHA")
present <- markers[markers %in% rownames(res_hyp)]
if (length(present) == 0)
  warning("No hypoxia marker genes found in results - check gene IDs / pipeline!")
print(as.data.frame(res_hyp[present, c("log2FoldChange", "padj")]))

cat("\nStep 3 complete. Five DE tables written to results/.\n")
