library(GEOquery)

dir.create("data", showWarnings = FALSE)
dir.create("metadata", showWarnings = FALSE)

# raw counts (one NCBI-generated matrix: Entrez gene IDs x samples)
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE101521&format=file&file=GSE101521_raw_counts_GRCh38.p13_NCBI.tsv.gz"
gz <- "data/GSE101521_counts.tsv.gz"
if (!file.exists(gz)) download.file(url, gz, mode = "wb")
counts <- as.matrix(read.delim(gzfile(gz), row.names = 1, check.names = FALSE))

# sample metadata; keep the HiSeq platform only (one sequencing batch)
gse <- getGEO("GSE101521", GSEMatrix = TRUE)
pd <- pData(gse[[grep("GPL16791", names(gse))]])

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

counts <- counts[, rownames(coldata)]

saveRDS(counts, "data/counts.rds")
saveRDS(coldata, "data/coldata.rds")
write.csv(data.frame(sample = rownames(coldata), coldata), "metadata/sample_sheet.csv", row.names = FALSE)

cat("genes x samples:", dim(counts)[1], "x", dim(counts)[2], "\n")
print(table(coldata$condition))
