#!/usr/bin/env nextflow

process SANKEY{
    /*
        This process takes outputs from previous processed
        and plots then in a sankey diagram per input type (read/contig)

        Inputs:
            - meta: Sample metadata
            - tobacco_mapping: path to mapping stats (reads to tobacco genome)
            - read2contigs: path to mapping stats (reads to contig fasta)
            - kraken2_read_summary_pluspf23: path to kraken2 output summary (read + pluspf23)
            - kraken2_read_summary_pluspf25: path to kraken2 output summary (read + pluspf25)
            - kraken2_read_summary_viper: path to kraken2 output summary (read + viper)
            - kraken2_contig_summary_pluspf23: path to kraken2 output summary (contig + pluspf23)
            - kraken2_contig_summary_pluspf25: path to kraken2 output summary (contig + pluspf25)
            - kraken2_contig_summary_viper: path to kraken2 output summary (contig + viper)

        Outputs:
            - read_sankey_output: tuple of meta with read sankey (html and json) and report text
            - contig_sankey_output: tuple of meta with contig sankey (html and json) and report text

    */
    container 'community.wave.seqera.io/library/mappy_pip_pandas_plotly:a119ac72e0c1a619'
    tag "${meta.id}"
    cpus 1
    memory '2GB'
    publishDir "${params.outdir}/${meta.id}/sankey", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(tobacco_mapping)
    tuple val(meta), path(read2contigs)
    tuple val(meta), path(kraken2_read_summary_pluspf23)
    tuple val(meta), path(kraken2_read_summary_pluspf25)
    tuple val(meta), path(kraken2_read_summary_viper)
    tuple val(meta), path(kraken2_contig_summary_pluspf23)
    tuple val(meta), path(kraken2_contig_summary_pluspf25)
    tuple val(meta), path(kraken2_contig_summary_viper)

    output:
    tuple val(meta), path("${meta.id}_read_sankey.html"), path("${meta.id}_read_sankey.json"), path("${meta.id}_read_sankey_text.txt"), emit: read_sankey_output
    tuple val(meta), path("${meta.id}_contig_sankey.html"), path("${meta.id}_contig_sankey.json"), path("${meta.id}_contig_sankey_text.txt"), emit: contig_sankey_output

    script:
    """
    unclassankey.py --sample ${meta.id} \
    --rank ${params.collapse_rank} \
    --tobacco_mapping ${tobacco_mapping} \
    --contigs ${read2contigs} \
    --kraken2_read_list ${kraken2_read_summary_pluspf23} ${kraken2_read_summary_pluspf25} ${kraken2_read_summary_viper} \
    --kraken2_contig_list ${kraken2_contig_summary_pluspf23} ${kraken2_contig_summary_pluspf25} ${kraken2_contig_summary_viper}
    """
}
