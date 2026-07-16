#!/usr/bin/env nextflow

process FLYE_OUTPUT_SUMMARY{
    /*
        This process takes outputs from flye, finds
        which reads are in whcih contigs and plots contig stats

        Inputs:
            - meta: Sample metadata
            - reads: Path to FASTQ file
            - assembly_fasta: Path to assembly.fasta output from flye
            - assembly_info: Path to assembly_info.txt output from flye

        Outputs:
            - read2contig_map: tuple of meta and which contig each read mapped to
            - contig_results: Tuple of meta and report outputs
    */
    container 'community.wave.seqera.io/library/mappy_pip_pandas_plotly:a119ac72e0c1a619'
    cpus 4
    memory '8GB'
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.id}/flye", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(gfa), path(assembly_fasta), path(assembly_info)

    output:
    tuple val(meta), path("${meta.id}_read2contig_mapping_stats.csv"), emit : read2contig_map
    tuple val(meta), path("${meta.id}_de_novo_assembly_stats_text.txt"), path("${meta.id}_de_novo_assembly_stats.html"), emit : contig_results

    script:
    """
    assembly_stats.py --sample ${meta.id} \
    --reads ${reads} \
    --assembly_fasta ${assembly_fasta} \
    --assembly_info ${assembly_info}
    """
}
