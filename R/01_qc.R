library(DESeq2)
library(pheatmap)

counts <- readRDS("data/counts.rds")
coldata <- readRDS("data/coldata.rds")
dir.create("figures", showWarnings = FALSE)

dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)

# drop genes with very low counts
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]

png("figures/libsize.png", width = 1000, height = 500)
barplot(colSums(counts(dds)), las = 2, cex.names = 0.5, main = "Library size")
dev.off()

vsd <- vst(dds)

png("figures/pca.png", width = 700, height = 600)
print(plotPCA(vsd, intgroup = "condition"))
dev.off()

d <- as.matrix(dist(t(assay(vsd))))
png("figures/sample_dist.png", width = 800, height = 700)
pheatmap(d, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = as.data.frame(colData(dds)["condition"]))
dev.off()

# flag and drop the clear PCA outlier (extreme PC1)
pc1 <- prcomp(t(assay(vsd)))$x[, 1]
outlier <- names(which.max(pc1))
cat("dropping outlier:", outlier, "\n")
dds <- dds[, colnames(dds) != outlier]

# PCA again without the outlier
png("figures/pca_clean.png", width = 700, height = 600)
print(plotPCA(vst(dds), intgroup = "condition"))
dev.off()

saveRDS(dds, "data/dds.rds")
cat("genes:", nrow(dds), " samples:", ncol(dds), "\n")
