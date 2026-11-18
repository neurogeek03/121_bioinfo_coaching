#!/usr/bin/env Rscript
library(GEOquery)                # parses GEO series-matrix files

# Both raw inputs are committed in the repo and staged into this work dir by
# Nextflow (see main.nf). We read them directly -- no network download, because
# the NCBI programmatic download is flaky. Filenames are bare: Nextflow stages
# the inputs here, and a direct `bin/00_download_data.R` run from the repo root
# finds the committed copies too.

## 1. Raw counts --------------------------------------------------------------
# NCBI-reprocessed count matrix: integer counts, Entrez gene IDs x samples.
# This file holds ALL 86 samples (both sequencing platforms) -- we subset later.
counts <- as.matrix(read.delim(gzfile("GSE101521_raw_counts_GRCh38.p13_NCBI.tsv.gz"),
  row.names = 1,            # 1st column (Entrez ID) -> row names
  check.names = FALSE))     # keep sample IDs verbatim (don't mangle to valid R names)

## 2. Sample metadata ---------------------------------------------------------
# GPL16791 (HiSeq) series matrix = 59 samples, one batch.
pd <- pData(getGEO(filename = "GSE101521-GPL16791_series_matrix.txt.gz", getGPL = FALSE))  # parse the LOCAL file; getGPL=FALSE skips the flaky network fetch of the platform annotation (we only need pData)

# ---------------------------------------------------------
# pause and inspect metadata
# ---------------------------------------------------------

# Pull out the fields we need. GEO stores everything as strings -> coerce numerics (change data type)
coldata <- data.frame(
  row.names = rownames(pd),                  # GSM sample IDs
  diagnosis = pd[["diagnosis:ch1"]],         # raw label, e.g. "major depression" (kept so the mapping is visible)
  age = as.numeric(pd[["age (yrs):ch1"]]),
  sex = pd[["Sex:ch1"]],
  rin = as.numeric(pd[["rin:ch1"]]),         # RNA integrity number
  pmi = as.numeric(pd[["pmi:ch1"]])          # post-mortem interval
)

# Collapse to a two-group contrast. CON first = DESeq2 reference level.
coldata$condition <- factor(
  ifelse(grepl("control", coldata$diagnosis, ignore.case = TRUE), "CON", "MDD"),
  levels = c("CON", "MDD"))
coldata$sex <- factor(coldata$sex)           # REQUIRED: the sex+condition model errors if sex is character

## 3. Align + save ------------------------------------------------------------
# One line that does two jobs: drop the non-HiSeq columns AND reorder the
# matrix to match coldata's rows. DESeq2 assumes counts cols == coldata rows.
counts <- counts[, rownames(coldata)]

saveRDS(counts, "counts.rds")                # binary handoff to 01_qc.R
saveRDS(coldata, "coldata.rds")
write.csv(data.frame(sample = rownames(coldata), coldata), "sample_sheet.csv", row.names = FALSE)  # human-readable copy

cat("genes x samples:", dim(counts)[1], "x", dim(counts)[2], "\n")   # sanity -> 39376 x 59
print(table(coldata$condition))              # -> CON 29, MDD 30
