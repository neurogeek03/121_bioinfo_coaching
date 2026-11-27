#!/usr/bin/env Rscript
library(DESeq2)
library(org.Hs.eg.db)            # human Entrez <-> gene symbol mapping

dds <- readRDS("dds.rds")        # filtered, outlier-dropped object from 01_qc.R

# The design formula is the pipeline's --design parameter. It IS the question we
# ask: "~ condition" = MDD vs CON; "~ sex + condition" = the same contrast, with
# sex held constant. Nothing else in the script changes.
design_str <- Sys.getenv("DESIGN", "~ condition")
design(dds) <- as.formula(design_str)
cat("design:", design_str, "\n")

# Columns we read in the results table, in the order we read them.
cols <- c("symbol", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")

## 1. Fit the model -----------------------------------------------------------
# One call does: size factors -> per-gene dispersion -> negative-binomial GLM.
dds <- DESeq(dds)

# lfcShrink shrinks noisy (low-count / high-variance) fold-changes toward 0 so
# the ranking is trustworthy. coef name comes from the CON reference level set in 00.
res <- lfcShrink(dds, coef = "condition_MDD_vs_CON", type = "apeglm")

## 2. Readable gene names -----------------------------------------------------
# Counts are keyed by Entrez IDs (stable; symbols get renamed over time). Add the
# symbol column so the results table is readable without a lookup.
res$symbol <- mapIds(org.Hs.eg.db, rownames(res), "SYMBOL", "ENTREZID")

## 3. Diagnostics -------------------------------------------------------------
# MA plot: fold-change (y) vs mean expression (x). Shrinkage visible at low counts.
# alpha = 0.05 so the coloured points are the same 45 genes we call significant
# everywhere else (plotMA's own default is 0.1, which would not match).
n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
png("ma.png", width = 760, height = 620)
plotMA(res, ylim = c(-3, 3), alpha = 0.05,
       main = "MA plot - effect size vs how much the gene is expressed",
       xlab = "A = Add  ->  average expression (mean of normalized counts)",
       ylab = "M = Minus  ->  log2(MDD) - log2(CON)")
legend("bottomright", bty = "n", cex = 0.95, pch = c(16, 16, 2),
       col = c("blue", "gray60", "gray40"),
       legend = c(sprintf("padj < 0.05  (%d genes)", n_sig),
                  "not significant",
                  "point clipped by the y-axis limit"))
dev.off()

# p-value histogram: should be ~uniform with a spike near 0 when signal exists.
png("pval_hist.png", width = 700, height = 500)
hist(res$pvalue, breaks = 50, col = "grey", main = "p-value distribution", xlab = "p-value")
dev.off()

## 4. Results table -----------------------------------------------------------
# Rank by padj (BH-adjusted); never judge on raw p across ~23k tests.
res <- res[order(res$padj), ]
write.csv(as.data.frame(res)[, cols], "dge_MDD_vs_CON.csv")
saveRDS(res, "res.rds")          # handoff to 03_interpret.R
cat("DE genes (", design_str, ") at padj < 0.05:",
    sum(res$padj < 0.05, na.rm = TRUE), "\n")   # ~ condition -> 45; ~ sex + condition -> 16

# Inspect what we just wrote: the top of the table, as we will read it together.
cat("\ntop 10 by padj:\n")
print(head(as.data.frame(res)[, cols], 10))
