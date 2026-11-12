#!/usr/bin/env nextflow

process WRITE_CONFIG{
    /*
        This process takes the template config file and adds the climb-id and dates.

        Inputs:
            - meta: Sample metadata
            - template: Path to template config file

        Outputs:
            - config: Tuple of meta edited config file
    */
    container 'community.wave.seqera.io/library/datetime_pip_argparse:59c77cbde267cdf7'
    tag "${meta.id}"
    cpus 2
    memory '2GB'
    publishDir "${params.outdir}/${meta.id}/report", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads), path(template)

    output:
    tuple val(meta), path("${meta.id}_config.json"), emit : config

    script:
    """
    write_config.py --sample ${meta.id} \
    --config $template \
    """
}