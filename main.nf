#!/usr/bin/env nextflow

// MDD vs control bulk RNA-seq: counts -> QC -> DGE -> interpretation.
// Each process wraps one R script in bin/. One teaching session = one step.

process DOWNLOAD {
    publishDir "${params.outdir}/data", mode: 'copy'

    input:
    path counts_gz          // raw counts, committed in the repo (NCBI download is flaky)
    path series_matrix_gz   // GPL16791 series matrix (sample metadata), committed in the repo

    output:
    path "counts.rds",  emit: counts
    path "coldata.rds", emit: coldata
    path "sample_sheet.csv"

    script:
    "00_load_data.R"
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
    // Raw inputs live in the repo; Nextflow stages them into the DOWNLOAD work dir.
    counts_gz        = Channel.fromPath("${projectDir}/data/GSE101521_raw_counts_GRCh38.p13_NCBI.tsv.gz")
    series_matrix_gz = Channel.fromPath("${projectDir}/data/GSE101521-GPL16791_series_matrix.txt.gz")

    DOWNLOAD(counts_gz, series_matrix_gz)
    QC(DOWNLOAD.out.counts, DOWNLOAD.out.coldata)
    DGE(QC.out.dds)
    INTERPRET(DGE.out.res, QC.out.dds)
}
