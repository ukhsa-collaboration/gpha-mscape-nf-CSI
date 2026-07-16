#!/usr/bin/env nextflow

process KRAKEN2_RESULT_SUMMARY {
    /*
        This process formats kraken stdout

        Inputs:
            - meta: Sample metadata
            - kraken2_stdout: Input file
            - database_name: val of database name
            - sampletype: val of reads/contigs


        Outputs:
            - kraken2_output_summary: tuple of meta and path to collapsed to family
            - kraken2_counts: tuple of meta and path to kraken2 results counted

    */
    executor 'local'
    conda 'assets/environment.yml'
    cpus 2
    memory '2GB'
    tag "${meta.id}"
    // errorStrategy 'ignore'
    publishDir "${params.outdir}/${meta.id}/kraken", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(database_name), val(sampletype), path(kraken2_stdout)

    output:
    tuple val(meta), path("${meta.id}_${sampletype}_kraken2_${database_name}_${params.collapse_rank}.csv"), emit: kraken2_output_summary
    tuple val(meta), path("${meta.id}_${sampletype}_kraken2_${database_name}_summary.csv"), emit: kraken2_counts

    script:
    """
    kraken_summary.py --sample ${meta.id} \
    --input_type ${sampletype} \
    --kraken_results ${kraken2_stdout} \
    --rank ${params.collapse_rank} \
    --database ${database_name}
    """
}
