#!/usr/bin/env nextflow

// MDD vs control bulk RNA-seq: counts -> QC -> DGE -> interpretation.
// Each process wraps one R script in bin/. One teaching session = one step.

process DOWNLOAD {
    publishDir "${params.outdir}/data", mode: 'copy'

    output:
    path "counts.rds",  emit: counts
    path "coldata.rds", emit: coldata
    path "sample_sheet.csv"

    script:
    "00_download_data.R"
}

process QC {
    publishDir "${params.outdir}/data",    mode: 'copy', pattern: "*.rds"
    publishDir "${params.outdir}/figures", mode: 'copy', pattern: "*.png"

    input:
    path "counts.rds"
    path "coldata.rds"

    output:
    path "dds.rds", emit: dds
    path "*.png"

    script:
    "01_qc.R"
}

process DGE {
    publishDir "${params.outdir}/data",    mode: 'copy', pattern: "res.rds"
    publishDir "${params.outdir}/figures", mode: 'copy', pattern: "*.png"
    publishDir "${params.outdir}/tables",  mode: 'copy', pattern: "*.csv"

    input:
    path "dds.rds"

    output:
    path "res.rds", emit: res
    path "*.csv"
    path "*.png"

    script:
    "02_dge.R"
}

process INTERPRET {
    publishDir "${params.outdir}/figures", mode: 'copy', pattern: "*.png"
    publishDir "${params.outdir}/tables",  mode: 'copy', pattern: "*.csv"

    input:
    path "res.rds"
    path "dds.rds"

    output:
    path "*.png"
    path "*.csv"

    script:
    "03_interpret.R"
}

workflow {
    DOWNLOAD()
    QC(DOWNLOAD.out.counts, DOWNLOAD.out.coldata)
    DGE(QC.out.dds)
    INTERPRET(DGE.out.res, QC.out.dds)
}
