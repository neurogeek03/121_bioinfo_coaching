#!/usr/bin/env Rscript           # run standalone; Nextflow calls it by name (bin/ is on PATH)
library(DESeq2)                  # count modelling + VST
library(pheatmap)               # the sample-distance heatmap

## 1. Build the DESeq2 object -------------------------------------------------
# Bundles counts + sample sheet + the model formula into one object. The files
# counts.rds / coldata.rds are staged in by Nextflow from the DOWNLOAD step.
counts <- readRDS("counts.rds")
coldata <- readRDS("coldata.rds")
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)

## 2. Filter low-count genes --------------------------------------------------
# Keep a gene only if >=10 samples have >=10 counts. Removes noise and shrinks
# the multiple-testing burden later (~39k -> ~23k genes).
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]

## 3. QC plots ----------------------------------------------------------------
# Sequencing depth per sample (wildly uneven depth is a red flag).
png("libsize.png", width = 1000, height = 500)
barplot(colSums(counts(dds)), las = 2, cex.names = 0.5, main = "Library size")
dev.off()

# VST = variance-stabilizing transform: puts samples on a comparable scale so
# plots aren't dominated by high-count genes. For visualization ONLY, not testing.
vsd <- vst(dds)

# PCA: do samples cluster by group? (For brain MDD, expect heavy overlap.)
png("pca.png", width = 700, height = 600)
print(plotPCA(vsd, intgroup = "condition"))
dev.off()

# Sample-to-sample distances: another view of outliers / batch structure.
d <- as.matrix(dist(t(assay(vsd))))
png("sample_dist.png", width = 800, height = 700)
pheatmap(d, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = as.data.frame(colData(dds)["condition"]))
dev.off()

## 4. Drop the PCA outlier ----------------------------------------------------
# One sample sits far out on PC1. Identify it, drop it, re-plot. Document the
# choice; don't rabbit-hole.
pc1 <- prcomp(t(assay(vsd)))$x[, 1]
outlier <- names(which.max(pc1))
cat("dropping outlier:", outlier, "\n")
dds <- dds[, colnames(dds) != outlier]

png("pca_clean.png", width = 700, height = 600)
print(plotPCA(vst(dds), intgroup = "condition"))
dev.off()

## 5. Save --------------------------------------------------------------------
saveRDS(dds, "dds.rds")          # binary handoff to 02_dge.R
cat("genes:", nrow(dds), " samples:", ncol(dds), "\n")   # sanity -> 23022 x 49
