#!/usr/bin/env nextflow

process RUN_KRAKEN2 {
    /*
        This process runs kraken on an input file.

        Inputs:
            - meta: Sample metadata
            - fasta/fastq: Input file

        Outputs:
            - kraken_outputs: Tuple of meta and kraken outputs
    */
    container 'quay.io/biocontainers/kraken2:2.17.1--pl5321h077b44d_0'
    cpus 8
    memory '10GB'
    tag "${meta.id}"
    // errorStrategy 'ignore'
    publishDir "${params.outdir}/${meta.id}/kraken", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(fasta)
    path(database)
    val(sampletype)
    val(database_name)

    output:
    tuple val(meta), val(database_name), val(sampletype), path("${meta.id}.${sampletype}.${database_name}.kraken_out.txt"), emit: kraken2_outputs

    script:
    """
    kraken2 --db $database \
    --report ${meta.id}.${sampletype}.${database_name}.kraken_report.txt \
    --threads 8 --memory-mapping $fasta > ${meta.id}.${sampletype}.${database_name}.kraken_out.txt
    """
}
