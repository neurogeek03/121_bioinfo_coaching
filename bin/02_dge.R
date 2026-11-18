#!/usr/bin/env Rscript
library(DESeq2)

dds <- readRDS("dds.rds")        # filtered, outlier-dropped object from 01_qc.R

## 1. Fit the model -----------------------------------------------------------
# One call does: size factors -> per-gene dispersion -> negative-binomial GLM.
dds <- DESeq(dds)

# lfcShrink shrinks noisy (low-count / high-variance) fold-changes toward 0 so
# the ranking is trustworthy. coef name comes from the CON reference level set in 00.
res <- lfcShrink(dds, coef = "condition_MDD_vs_CON", type = "apeglm")

## 2. Diagnostics -------------------------------------------------------------
# MA plot: fold-change vs mean expression (shrinkage visible at low counts).
png("ma.png", width = 700, height = 600)
plotMA(res, ylim = c(-3, 3))
dev.off()

# p-value histogram: should be ~uniform with a spike near 0 when signal exists.
png("pval_hist.png", width = 700, height = 500)
hist(res$pvalue, breaks = 50, col = "grey", main = "p-value distribution", xlab = "p-value")
dev.off()

## 3. Results table -----------------------------------------------------------
# Rank by padj (BH-adjusted); never judge on raw p across ~23k tests.
res <- res[order(res$padj), ]
write.csv(as.data.frame(res), "dge_MDD_vs_CON.csv")
saveRDS(res, "res.rds")          # handoff to 03_interpret.R
cat("DE genes (~ condition) at padj < 0.05:", sum(res$padj < 0.05, na.rm = TRUE), "\n")   # -> 45

## 4. Covariate model ---------------------------------------------------------
# Adjust for sex: the formula IS the experimental design. Controlling for sex
# removes sex-confounded signal and costs power -> fewer hits (45 -> 16).
design(dds) <- ~ sex + condition
dds <- DESeq(dds)
res_sex <- lfcShrink(dds, coef = "condition_MDD_vs_CON", type = "apeglm")
res_sex <- res_sex[order(res_sex$padj), ]
write.csv(as.data.frame(res_sex), "dge_MDD_vs_CON_adjSex.csv")
cat("DE genes (~ sex + condition) at padj < 0.05:", sum(res_sex$padj < 0.05, na.rm = TRUE), "\n")   # -> 16
