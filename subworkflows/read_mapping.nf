#!/usr/bin/env nextflow

include { RUN_MAPPING } from '../modules/mapping'

workflow READ_MAPPING {
    take:
    ch_reads   // channel: [ val(meta), path(reads), path(refs)]

    main:
        RUN_MAPPING(ch_reads)

    emit:
    mappy_stats = RUN_MAPPING.out.mappy_stats  // channel: [ val(meta), path(mapping_stats) ]
    mappy_results = RUN_MAPPING.out.mappy_results   // channel: [ val(meta), path(mapping_stats_text) ]

}
