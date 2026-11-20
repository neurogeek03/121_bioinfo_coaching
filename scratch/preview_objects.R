#!/usr/bin/env Rscript
# Session 2 opener — preview the objects the DOWNLOAD step built, BEFORE QC
# feeds them to DESeq2. Nothing here changes the data; it just looks at it.
#
# Run from the repo root:   Rscript scratch/preview_objects.R
# or open in RStudio and step through line by line (nicer for teaching).

counts  <- readRDS("results/data/counts.rds")   # integer matrix: genes x samples
coldata <- readRDS("results/data/coldata.rds")  # the sample sheet (design / metadata)

## 1. What ARE these two objects? --------------------------------------------
class(counts)                 # "matrix" "array"   -- one grid, all one type
class(coldata)                # "data.frame"        -- columns of mixed types

## 2. The counts matrix -------------------------------------------------------
dim(counts)                   # 39376 x 59  (ALL gene biotypes, before filtering)
head(counts[, 1:6])           # first 6 genes x 6 samples: integer counts,
                              #   Entrez IDs as rownames, GSM IDs as colnames

## 3. The sample sheet --------------------------------------------------------
str(coldata)                  # factors stayed factors, numerics stayed numeric
head(coldata)                 # first 6 samples and all their variables
table(coldata$condition)      # CON 29 / MDD 30

## 4. Do the two line up? (why the align line in step 00 mattered) ------------
all(colnames(counts) == rownames(coldata))   # TRUE: counts cols == coldata rows

## 5. Missing values (the "NAs introduced by coercion" warning from step 00) --
colSums(is.na(coldata))       # rin 2, pmi 1 -- everything else clean
head(coldata[!complete.cases(coldata), ])    # the exact samples with a gap
