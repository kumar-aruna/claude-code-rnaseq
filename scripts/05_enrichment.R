# ============================================================
# Step 5: Pathway / enrichment analysis (via enrichR / Enrichr)
# Goal: turn gene LISTS into BIOLOGY. For a set of changed genes,
#       find which known pathways / gene sets are over-represented.
# Tool: enrichR queries the Enrichr web service (needs internet,
#       no compilation needed). Databases used:
#         - MSigDB_Hallmark_2020   (classic cancer/biology signatures)
#         - GO_Biological_Process_2021
#         - KEGG_2021_Human        (canonical pathways)
# Output: enrichment tables (CSV) + barplots (PNG) in results/
# ============================================================

suppressMessages({
  library(enrichR)
  library(ggplot2)
})

setEnrichrSite("Enrichr")                         # human gene sets
DBS <- c("MSigDB_Hallmark_2020",
         "GO_Biological_Process_2021",
         "KEGG_2021_Human")

# --- Helper: load a DE table, return up/down/all significant genes ---
load_sets <- function(csv) {
  if (!file.exists(csv))
    stop("Missing ", csv, " - run scripts/03_deseq.R first.")
  d <- read.csv(csv, row.names = 1)
  d <- d[!is.na(d$padj), ]
  list(up   = rownames(d[d$padj < 0.05 & d$log2FoldChange >  1, ]),
       down = rownames(d[d$padj < 0.05 & d$log2FoldChange < -1, ]),
       all  = rownames(d[d$padj < 0.05 & abs(d$log2FoldChange) > 1, ]))
}

# --- Helper: run Enrichr on a gene set, save tables + a barplot ------
enrich_one <- function(genes, label, prefix) {
  cat(sprintf("\n%s  (%d genes)\n", label, length(genes)))
  if (length(genes) < 5) { cat("  too few genes, skipping\n"); return(invisible()) }

  res <- enrichr(genes, DBS)                       # list: one table per DB

  for (db in names(res)) {                         # save every DB's table
    tab <- res[[db]]
    if (nrow(tab) > 0)
      write.csv(tab, sprintf("%s_%s.csv", prefix, db), row.names = FALSE)
  }

  # Barplot of the top 12 Hallmark gene sets (most interpretable view)
  hm <- res[["MSigDB_Hallmark_2020"]]
  hm <- hm[order(hm$Adjusted.P.value), ]
  hm <- head(hm[hm$Adjusted.P.value < 0.05, ], 12)
  if (nrow(hm) > 0) {
    hm$Term <- factor(hm$Term, levels = rev(hm$Term))
    p <- ggplot(hm, aes(-log10(Adjusted.P.value), Term)) +
      geom_col(fill = "firebrick") +
      labs(title = paste("Hallmark enrichment:", label),
           x = "-log10(adjusted p)", y = NULL) +
      theme_bw()
    ggsave(sprintf("%s_hallmark.png", prefix), p, width = 9, height = 5, dpi = 150)
    # rows are ordered by ascending adjusted p, so row 1 is the strongest term
    cat("  top Hallmark term:", as.character(hm$Term[1]),
        "| saved", sprintf("%s_hallmark.png\n", prefix))
  } else {
    cat("  no significant Hallmark terms\n")
  }
  invisible(res)
}

# 1. Hypoxia UP genes -> expect HALLMARK_HYPOXIA / glycolysis on top.
s1 <- load_sets("results/DE_1_hypoxia_vs_normoxia_in_control.csv")
enrich_one(s1$up, "Hypoxia UP genes (Control)", "results/enrich_hypoxia_up")

# 2. ITPR3 x hypoxia interaction genes.
s4 <- load_sets("results/DE_4_ITPR3_hypoxia_interaction.csv")
enrich_one(s4$all, "ITPR3 x hypoxia interaction", "results/enrich_ITPR3_interaction")

# 3. RELB x hypoxia interaction genes.
s5 <- load_sets("results/DE_5_RELB_hypoxia_interaction.csv")
enrich_one(s5$all, "RELB x hypoxia interaction", "results/enrich_RELB_interaction")

cat("\nStep 5 complete.\n")
