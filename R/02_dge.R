library(DESeq2)

dds <- readRDS("data/dds.rds")
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

dds <- DESeq(dds)
res <- lfcShrink(dds, coef = "condition_MDD_vs_CON", type = "apeglm")

png("figures/ma.png", width = 700, height = 600)
plotMA(res, ylim = c(-3, 3))
dev.off()

png("figures/pval_hist.png", width = 700, height = 500)
hist(res$pvalue, breaks = 50, col = "grey", main = "p-value distribution", xlab = "p-value")
dev.off()

res <- res[order(res$padj), ]
write.csv(as.data.frame(res), "results/dge_MDD_vs_CON.csv")
saveRDS(res, "data/res.rds")

cat("DE genes (~ condition) at padj < 0.05:", sum(res$padj < 0.05, na.rm = TRUE), "\n")

# covariate model: adjust for sex
design(dds) <- ~ sex + condition
dds <- DESeq(dds)
res_sex <- lfcShrink(dds, coef = "condition_MDD_vs_CON", type = "apeglm")
res_sex <- res_sex[order(res_sex$padj), ]
write.csv(as.data.frame(res_sex), "results/dge_MDD_vs_CON_adjSex.csv")

cat("DE genes (~ sex + condition) at padj < 0.05:", sum(res_sex$padj < 0.05, na.rm = TRUE), "\n")
