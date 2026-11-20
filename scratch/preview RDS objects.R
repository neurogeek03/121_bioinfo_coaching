#!/usr/bin/env Rscript
# Session 2 opener — LOOK at the objects the DOWNLOAD step built, before QC
# feeds them to DESeq2. Nothing here changes anything; it just inspects them.
#
# Run from the repo root:   Rscript "scratch/preview RDS objects.R"
# or open in RStudio / VS Code and step through it line by line (nicer for teaching).
# Companion script: "preview DDS object.R" inspects the dds the QC step builds.

counts  <- readRDS("counts.rds")   # integer matrix: genes x samples
coldata <- readRDS("coldata.rds")  # the sample sheet (design / metadata)

## ---- The counts matrix: what does it look like? ----------------------------
head(counts[, 1:6])           # <- first 6 genes x 6 samples, so it fits the screen
                              #    integer counts; Entrez IDs = rownames, GSM = colnames
class(counts)                 # "matrix" "array"  -- one grid, all one type
dim(counts)                   # 39376 x 59        -- ALL gene biotypes, before filtering

## ---- The sample sheet: what does it look like? -----------------------------
head(coldata)                 # <- first 6 samples and all their variables
class(coldata)                # "data.frame"      -- columns of mixed types
str(coldata)                  # factors stayed factors, numerics stayed numeric
table(coldata$condition)      # CON 29 / MDD 30

## ---- Do the two line up? (why the align line in step 00 mattered) ----------
all(colnames(counts) == rownames(coldata))   # TRUE: counts cols == coldata rows

## ---- Missing values (the "NAs introduced by coercion" warning from step 00) -
colSums(is.na(coldata))       # rin 2, pmi 1 -- everything else clean
head(coldata[!complete.cases(coldata), ])    # the exact samples with a gap
