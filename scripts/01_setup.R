# ============================================================
# Step 1: Setup & sanity checks
# Project: GSE197576 — ITPR3 / RELB knockout in SW480, normoxia vs hypoxia
# Goal of this script:
#   1. Load the raw gene count matrix
#   2. Build a sample metadata table (which sample = which condition)
#   3. Sanity-check that everything lines up
#   4. Build a DESeq2 object, ready for the next step
# ============================================================

# --- Load the packages we need ---
library(DESeq2)

# --- 1. Load the raw count matrix ------------------------------------
# The file has genes as rows and 12 samples as columns.
# row.names = 1  -> use the first column ("gene") as row names.
# check.names = FALSE -> keep sample names exactly as written.
# All paths are relative to the project root, so run scripts from there:
#   Rscript scripts/01_setup.R   (not from inside scripts/)
data_file <- "data/GSE197576_raw_gene_counts_matrix.tsv"
if (!file.exists(data_file))
  stop("Cannot find ", data_file,
       " - run this from the project root (cd into rnaseq-analysis first).")
dir.create("results", showWarnings = FALSE)   # ensure output dir exists

counts <- read.delim(data_file, row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)

cat("Count matrix loaded:\n")
cat("  genes (rows):   ", nrow(counts), "\n")
cat("  samples (cols): ", ncol(counts), "\n\n")

# --- 2. Build the sample metadata table ------------------------------
# We DECODE the condition from each column name instead of typing it by
# hand, so there is no risk of a mismatch. The column names look like:
#   01_SW_sgCTRL_Norm  ->  genotype = CTRL,  oxygen = Norm
sample_names <- colnames(counts)

# genotype: pull out CTRL / ITPR3 / RELB from the "sg..." part
genotype <- ifelse(grepl("sgCTRL",  sample_names), "Control",
             ifelse(grepl("sgITPR3", sample_names), "ITPR3_KO",
             ifelse(grepl("sgRELB",  sample_names), "RELB_KO", NA)))

# oxygen: Norm -> Normoxia, Hyp -> Hypoxia
oxygen <- ifelse(grepl("_Norm", sample_names), "Normoxia",
           ifelse(grepl("_Hyp",  sample_names), "Hypoxia", NA))

# Make them factors and set the REFERENCE level explicitly.
# DESeq2 compares everything *against* the reference, so Control and
# Normoxia are the natural baselines.
coldata <- data.frame(
  row.names = sample_names,
  genotype  = factor(genotype, levels = c("Control", "ITPR3_KO", "RELB_KO")),
  oxygen    = factor(oxygen,   levels = c("Normoxia", "Hypoxia"))
)

cat("Sample metadata table:\n")
print(coldata)
cat("\n")

# --- 3. Sanity checks ------------------------------------------------
# These stops() will halt the script loudly if anything is wrong,
# which is exactly what we want before trusting any result.
stopifnot(
  "No NA allowed in genotype" = !any(is.na(coldata$genotype)),
  "No NA allowed in oxygen"   = !any(is.na(coldata$oxygen)),
  # The metadata rows MUST be in the same order as the matrix columns:
  "coldata rows must match count columns" =
    all(rownames(coldata) == colnames(counts))
)
cat("Sanity checks passed: metadata matches the 12 samples.\n\n")

cat("Design summary (samples per group):\n")
print(table(coldata$genotype, coldata$oxygen))
cat("\n")

# --- 4. Build the DESeq2 object --------------------------------------
# The design formula ~ genotype * oxygen models:
#   - the effect of genotype
#   - the effect of oxygen
#   - the INTERACTION (does a knockout change the hypoxia response?)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ genotype * oxygen
)

# Light pre-filter: drop genes with almost no reads across all samples.
# This speeds things up and removes noise. (Full QC comes in Step 2.)
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat("Genes kept after light filter (rowSums >= 10):", nrow(dds), "\n\n")

# Save the object so the next step can pick up where we left off.
saveRDS(dds, "results/dds_step1.rds")
cat("Saved DESeq2 object to results/dds_step1.rds\n")
cat("Step 1 complete.\n")
