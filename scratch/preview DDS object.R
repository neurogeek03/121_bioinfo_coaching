#!/usr/bin/env Rscript
# Session 2 — LOOK at the dds (DESeqDataSet) the QC step built. The RDS preview
# script showed two loose tables (counts + coldata); this shows how DESeq2
# bundles them into ONE object that keeps them aligned. Nothing here changes
# anything; it just inspects.
#
# Run from the repo root:   Rscript "scratch/preview DDS object.R"
# or open in RStudio / VS Code and step through it line by line (nicer for teaching).
# Companion script: "preview RDS objects.R" previews the raw counts + sample sheet.

library(DESeq2)                        # needed so the show method + accessors work

dds <- readRDS("dds.rds")              # the object QC saved (post-filter, post-outlier-drop)

## ---- Print it: the show method -------------------------------------------
dds                                    # one screen: class, dim, assays, rownames/colnames, colData names

## ---- What KIND of object is it? ------------------------------------------
class(dds)                             # "DESeqDataSet"  -- not a data.frame
isS4(dds)                              # TRUE            -- S4: formal named slots
dim(dds)                               # 23292 x 58      -- genes x samples (after filter + outlier)

## ---- The three coordinated tables inside ---------------------------------
counts(dds)[1:3, 1:3]                  # the integer count matrix (genes x samples)
colData(dds)                           # the sample sheet -- the COLUMNS of the counts
dds$condition                          # pull one colData column directly
rowData(dds)                           # per-gene table   -- the ROWS of the counts

## ---- The model (the question, not data) ----------------------------------
design(dds)                            # ~ condition  -- "explain each gene's counts by condition"
assayNames(dds)                        # which matrices live inside (grows after vst/DESeq)
sizeFactors(dds)                       # NULL until DESeq() runs -- nothing fitted yet

## ---- Why bundle? the coordinates stay aligned ----------------------------
# One index subsets counts AND colData together -- nothing drifts out of sync.
mdd <- dds[, dds$condition == "MDD"]   # keep MDD samples only
dim(mdd)                               # fewer columns; counts cols == colData rows, still matched
