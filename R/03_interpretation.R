library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

res <- readRDS("data/res.rds")
dds <- readRDS("data/dds.rds")
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# Entrez ID -> gene symbol
res$symbol <- mapIds(org.Hs.eg.db, rownames(res), "SYMBOL", "ENTREZID")

png("figures/volcano.png", width = 900, height = 800)
print(EnhancedVolcano(res, lab = res$symbol,
                      x = "log2FoldChange", y = "padj",
                      pCutoff = 0.05, FCcutoff = 1))
dev.off()

top <- head(rownames(res), 30)
vsd <- vst(dds)
mat <- assay(vsd)[top, ]
rownames(mat) <- res[top, "symbol"]
png("figures/heatmap_top.png", width = 700, height = 800)
pheatmap(mat, scale = "row", show_colnames = FALSE,
         annotation_col = as.data.frame(colData(dds)["condition"]))
dev.off()

# GO enrichment on significant genes
sig <- rownames(res)[which(res$padj < 0.05)]
ego <- enrichGO(sig, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                ont = "BP", universe = rownames(res), readable = TRUE)
write.csv(as.data.frame(ego), "results/go_BP.csv", row.names = FALSE)

cat("significant genes:", length(sig), " | enriched GO terms:", nrow(ego), "\n")
