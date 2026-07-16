#!/usr/bin/env nextflow

process RUN_MAPPING{
    /*
        This process takes reads from a sample and maps
        them to a reference sequence using mappy.

        Inputs:
            - meta: Sample metadata
            - taxon_reads: Path to FASTQ file
            - ref_fasta: Path to reference FASTA to map against

        Outputs:
            - mappy_stats: Tuple of meta and mappy stats per read
            - mappy_results: Tuple of meta and report outputs
    */
    container 'community.wave.seqera.io/library/mappy_pip_pandas_plotly:a119ac72e0c1a619'
    tag "${meta.id}"
    cpus 8
    memory '32GB'
    publishDir "${params.outdir}/${meta.id}/mapping", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads), path(refs)

    output:
    tuple val(meta), path("${meta.id}_mapping_stats.csv"), emit : mappy_stats
    tuple val(meta), path("${meta.id}_mapping_stats_text.txt"), emit: mappy_results // will need changing when you want pileups etc.

    script:
    """
    mapping.py --sample ${meta.id} \
    --fastq $reads \
    --references $refs
    """
}
