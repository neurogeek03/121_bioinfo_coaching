# Bulk RNA-seq coaching — MDD vs control (GSE101521)

Four-session intro to a bulk RNA-seq case–control analysis in R: QC → differential
expression → interpretation. Dataset: postmortem dorsolateral PFC, major depressive
disorder vs non-psychiatric controls (Pantazatos et al., *Mol Psychiatry* 2017,
[GSE101521](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101521), PMID 27528462).

## Environment
Open in GitHub Codespaces. The dev container (`.devcontainer/devcontainer.json`) uses the
Bioconductor image and installs the packages on first build — no local R setup needed.

## Run order
```
Rscript R/00_download_data.R   # download counts, build sample sheet (29 CON vs 30 MDD)
Rscript R/01_qc.R              # filtering, VST, PCA, sample-distance heatmap
Rscript R/02_dge.R             # DESeq2 model, shrinkage, results table, MA + p-value plots
Rscript R/03_interpretation.R  # volcano, top-gene heatmap, GO enrichment
```
Each script saves its objects under `data/` and writes plots to `figures/` and tables to
`results/`, so a later step can start from the previous one.

## Sessions
1. Setup + data into R (counts matrix, sample sheet, factors)
2. QC (library sizes, filtering, normalization, PCA, outliers)
3. Differential expression (design, DESeq2, shrinkage, reading results)
4. Interpretation (volcano, heatmap, GO enrichment, biology)

## Covariates available
`age`, `sex`, `rin`, `pmi` — for the model-design extension (`~ sex + condition`, etc.).
