#!/usr/bin/env nextflow

process TOBACCO_CONTAMINATION_CHECK{
    /*
        This process takes tobacco mapping results.
        If >70 reads mapped to genome - it is classed as
        tobacco contamination (True).

        Inputs:
            - meta: Sample metadata
            - mappy_text: Text summary of mapping results

        Outputs:
            - tobacco_contamination: txt with True (contaminated -> will end pipeline) or False (not contaminated -> will continue pipeline)
    */
    tag "${meta.id}"
    cpus 1
    memory '2GB'
    container 'community.wave.seqera.io/library/mappy_pip_pandas_plotly:a119ac72e0c1a619'

    input:
    tuple val(meta), path(mappy_text)
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path(fastq), path("tobacco_contamination_check.txt"), emit : tobacco_contamination_results

    script:
    """
    PCT=\$(grep -r \\% ${mappy_text} | cut -d ' ' -f 4 | grep -Eo '[[:digit:]]+([.][[:digit:]]+)?')

    if [[ -z "\$PCT" ]]; then
        PCT=0
    fi

    if [[ \$PCT > 70 ]]; then
        result="True"
        echo \$result > tobacco_contamination_check.txt
    else
        result="False"
        echo \$result > tobacco_contamination_check.txt
    fi
    """

}
