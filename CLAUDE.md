# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project status

Bulk RNA-seq analysis in **R + DESeq2** (v1.46.0). Analysis is built as numbered
scripts in `scripts/`, run in order from the project root. Each script saves an
intermediate object or figures so the next step can resume.

Progress: ✅ Step 1 setup, ✅ Step 2 QC. 🔜 Step 3 differential expression.

### How to run

```bash
Rscript scripts/01_setup.R   # load data, build metadata, make DESeqDataSet -> results/dds_step1.rds
Rscript scripts/02_qc.R      # VST + PCA + sample-distance heatmap -> results/qc_*.png
```

### Environment notes

- Install bio packages via `BiocManager::install(...)`, NOT `install.packages()`.
- `apeglm`/`ashr` are NOT installed (no gfortran on this machine), so use
  `lfcShrink(..., type = "normal")` for fold-change shrinkage — do not use the
  apeglm default.

## Directory layout

```
data/        # input data (count matrix + gzipped backup); treat as read-only
scripts/     # numbered analysis scripts (01_setup, 02_qc, ...)
results/     # generated objects (.rds), figures (.png), DE tables (.csv)
notebooks/   # exploratory work (empty for now)
```

## Dataset: GSE197576

The starting point is a pre-computed raw gene count matrix from GEO accession **GSE197576**, so the FASTQ → alignment → quantification steps are already done. Work begins at differential expression.

- File: `data/GSE197576_raw_gene_counts_matrix.tsv` (gzipped backup kept alongside).
- Shape: **43,809 genes × 12 samples**. First column is `gene`; remaining 12 are samples.
- Gene IDs are **HGNC symbols** (not Ensembl) — no ID conversion needed.
- Values are **raw integer counts** — correct input for DESeq2/edgeR/pydeseq2. Filter low-count genes before testing.
- Organism: *Homo sapiens*. Cell line: SW480 (colorectal). Platform: Illumina NextSeq 500. BioProject PRJNA811111.

## Experimental design (2 × 3 factorial, 2 replicates)

The study knocks out **ITPR3** or **RELB** vs. control, under **normoxia** vs. **hypoxia** (0.5% O₂, 72h). Sample column names encode genotype and oxygen condition; decode them as:

| Genotype | Normoxia | Hypoxia |
|----------|----------|---------|
| Control (`sgCTRL`) | `01_SW_sgCTRL_Norm`, `02_SW_sgCTRL_Norm` | `11_SW_sgCTRL_Hyp`, `12_SW_sgCTRL_Hyp` |
| ITPR3 KO (`sgITPR3_1`) | `03_SW_sgITPR3_1_Norm`, `04_SW_sgITPR3_1_Norm` | `13_SW_sgITPR3_1_Hyp`, `14_SW_sgITPR3_1_Hyp` |
| RELB KO (`sgRELB_3`) | `07_SW_sgRELB_3_Norm`, `08_SW_sgRELB_3_Norm` | `17_SW_sgRELB_3_Hyp`, `18_SW_sgRELB_3_Hyp` |

Notes that matter for modeling:
- Two factors (**genotype**, **oxygen**) with 2 biological replicates per cell — supports a `~ genotype * oxygen` interaction model.
- Column numbering is non-contiguous (01–18 with gaps) because samples from the full original study were subset; do not assume sequential IDs.
- Build the sample metadata table by parsing the column headers, and keep factor level ordering explicit (e.g. control and normoxia as reference levels).
