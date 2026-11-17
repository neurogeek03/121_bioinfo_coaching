#!/usr/bin/env Rscript
library(GEOquery)

# raw counts (one NCBI-generated matrix: Entrez gene IDs x samples)
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE101521&format=file&file=GSE101521_raw_counts_GRCh38.p13_NCBI.tsv.gz"
download.file(url, "counts.tsv.gz", mode = "wb")
counts <- as.matrix(read.delim(gzfile("counts.tsv.gz"), row.names = 1, check.names = FALSE))

# sample metadata; keep the HiSeq platform only (one sequencing batch)
mtx <- "GSE101521-GPL16791_series_matrix.txt.gz"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE101nnn/GSE101521/matrix/GSE101521-GPL16791_series_matrix.txt.gz", mtx, mode = "wb")
pd <- pData(getGEO(filename = mtx))

coldata <- data.frame(
  row.names = rownames(pd),
  diagnosis = pd[["diagnosis:ch1"]],
  age = as.numeric(pd[["age (yrs):ch1"]]),
  sex = pd[["Sex:ch1"]],
  rin = as.numeric(pd[["rin:ch1"]]),
  pmi = as.numeric(pd[["pmi:ch1"]])
)

# collapse to a simple two-group contrast: control vs MDD
coldata$condition <- ifelse(grepl("control", coldata$diagnosis, ignore.case = TRUE), "CON", "MDD")
coldata$condition <- factor(coldata$condition, levels = c("CON", "MDD"))
coldata$sex <- factor(coldata$sex)

counts <- counts[, rownames(coldata)]

saveRDS(counts, "counts.rds")
saveRDS(coldata, "coldata.rds")
write.csv(data.frame(sample = rownames(coldata), coldata), "sample_sheet.csv", row.names = FALSE)

cat("genes x samples:", dim(counts)[1], "x", dim(counts)[2], "\n")
print(table(coldata$condition))
