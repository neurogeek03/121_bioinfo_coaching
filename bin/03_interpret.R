#!/usr/bin/env Rscript
library(DESeq2)
library(EnhancedVolcano)        # volcano plot
library(pheatmap)               # top-gene heatmap
library(clusterProfiler)        # GO over-representation
library(org.Hs.eg.db)           # human Entrez <-> symbol / GO mapping

res <- readRDS("res.rds")        # results from 02_dge.R
dds <- readRDS("dds.rds")        # pre-DESeq object from 01_qc.R (for the log transform + colData)

## 1. Readable gene names -----------------------------------------------------
# Results are keyed by Entrez IDs (stable, from the upstream count pipeline);
# map to symbols so a human can read the plots.
res$symbol <- mapIds(org.Hs.eg.db, rownames(res), "SYMBOL", "ENTREZID")

## 2. Volcano ----------------------------------------------------------------
# Effect size (log2FC) vs significance (padj). Top-right/left = the interesting corner.
png("volcano.png", width = 900, height = 800)
print(EnhancedVolcano(res, lab = res$symbol,
                      x = "log2FoldChange", y = "padj",
                      pCutoff = 0.05, FCcutoff = 1))
dev.off()

## 3. Heatmap of the top genes ------------------------------------------------
# scale="row" shows each gene's relative pattern across samples (z-scored).
top <- head(rownames(res), 30)
ntd <- normTransform(dds)        # log2(normalized + 1), same transform as QC
mat <- assay(ntd)[top, ]
rownames(mat) <- res[top, "symbol"]
png("heatmap_top.png", width = 700, height = 800)
pheatmap(mat, scale = "row", show_colnames = FALSE,
         annotation_col = as.data.frame(colData(dds)["condition"]))
dev.off()

## 4. GO enrichment -----------------------------------------------------------
# What biology do the significant genes share? universe = expressed genes tested
# (NOT the whole genome) or enrichment is inflated.
sig <- rownames(res)[which(res$padj < 0.05)]
ego <- enrichGO(sig, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                ont = "BP", universe = rownames(res), readable = TRUE)
write.csv(as.data.frame(ego), "go_BP.csv", row.names = FALSE)
cat("significant genes:", length(sig), " | enriched GO terms:", nrow(ego), "\n")   # -> 45 | 458
