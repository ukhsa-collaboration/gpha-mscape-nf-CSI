#!/usr/bin/env nextflow

include { SANKEY } from '../modules/sankey'

workflow PLOT_CLASSIFICATION_RESULTS {
    take:
    ch_all_results

    main:
        SANKEY(ch_all_results)

}
