#!/usr/bin/env Rscript
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

res <- readRDS("res.rds")
dds <- readRDS("dds.rds")

# Entrez ID -> gene symbol
res$symbol <- mapIds(org.Hs.eg.db, rownames(res), "SYMBOL", "ENTREZID")

png("volcano.png", width = 900, height = 800)
print(EnhancedVolcano(res, lab = res$symbol,
                      x = "log2FoldChange", y = "padj",
                      pCutoff = 0.05, FCcutoff = 1))
dev.off()

top <- head(rownames(res), 30)
vsd <- vst(dds)
mat <- assay(vsd)[top, ]
rownames(mat) <- res[top, "symbol"]
png("heatmap_top.png", width = 700, height = 800)
pheatmap(mat, scale = "row", show_colnames = FALSE,
         annotation_col = as.data.frame(colData(dds)["condition"]))
dev.off()

# GO enrichment on significant genes
sig <- rownames(res)[which(res$padj < 0.05)]
ego <- enrichGO(sig, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                ont = "BP", universe = rownames(res), readable = TRUE)
write.csv(as.data.frame(ego), "go_BP.csv", row.names = FALSE)
cat("significant genes:", length(sig), " | enriched GO terms:", nrow(ego), "\n")
