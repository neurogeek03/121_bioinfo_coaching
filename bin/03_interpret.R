#!/usr/bin/env Rscript
library(DESeq2)
library(EnhancedVolcano)        # volcano plot
library(clusterProfiler)        # GO over-representation
library(enrichplot)             # dotplot of the enriched terms
library(ggplot2)                # ggtitle on the dotplot
library(org.Hs.eg.db)           # human Entrez <-> symbol / GO mapping

res <- readRDS("res.rds")        # results from 02_dge.R

# res already carries a `symbol` column, mapped from Entrez IDs in 02_dge.R.

## 1. Volcano ----------------------------------------------------------------
# Effect size (log2FC) vs significance (padj). Top-right/left = the interesting corner.
png("volcano.png", width = 900, height = 800)
print(EnhancedVolcano(res, lab = res$symbol,
                      x = "log2FoldChange", y = "padj",
                      pCutoff = 0.05, FCcutoff = 1))
dev.off()

## 2. GO enrichment -----------------------------------------------------------
# What biology do the significant genes share? universe = expressed genes tested
# (NOT the whole genome) or enrichment is inflated.
sig <- rownames(res)[which(res$padj < 0.05)]
ego <- enrichGO(sig, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                ont = "BP", universe = rownames(res), readable = TRUE)
write.csv(as.data.frame(ego), "go_BP.csv", row.names = FALSE)

# Dotplot of the top terms: y = GO term, x = GeneRatio (what share of our
# significant genes sit in that term), dot size = gene count, colour = p.adjust.
# This is the figure to read out loud -- the CSV is for looking things up.
png("go_dotplot.png", width = 900, height = 700, res = 110)
print(dotplot(ego, showCategory = 15) + ggtitle("Enriched GO biological processes"))
dev.off()

cat("significant genes:", length(sig), " | enriched GO terms:", nrow(ego), "\n")   # -> 45 | 458
