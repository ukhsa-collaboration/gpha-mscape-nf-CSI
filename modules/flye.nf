#!/usr/bin/env nextflow

process RUN_FLYE {
    /*
        This process takes a fastq file and generates a de novo assembly. Setting for
        flye are optimised for subsequent use of strainy to phase genomes.

        Inputs:
            - meta: Sample metadata
            - fastq: Input FASTQ file

        Outputs:
            - flye_outputs: Tuple of meta and flye output files
            - flye_logs
    */
    container 'quay.io/biocontainers/flye:2.9.6--py310h5850263_1'
    tag "${meta.id}"
    cpus 8
    memory '10GB'
    // errorStrategy 'ignore'
    publishDir "${params.outdir}/${meta.id}/flye", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("assembly_graph.gfa"), path("assembly.fasta"), path("assembly_info.txt"), emit: flye_outputs
    tuple val(meta), path("flye.log"), emit: flye_logs

    script:
    """
    flye --nano-raw $fastq -o . -t 8 --meta --no-alt-contigs -i 0
    """
}
