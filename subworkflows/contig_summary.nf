#!/usr/bin/env nextflow

include { FLYE_OUTPUT_SUMMARY } from '../modules/flye_summary'

workflow CONTIG_SUMMARY {
    take:
        ch_reads_pass   // channel: [ val(meta), path(reads)] - pass tobacco contamination check
        ch_contigs // channel: [ val(meta), val(gfa), val(assembly), val(info)]

    main:
        FLYE_OUTPUT_SUMMARY(ch_reads_pass,
                            ch_contigs)

    emit:
    read2contigs = FLYE_OUTPUT_SUMMARY.out.read2contig_map  // channel: [ val(meta), path(csv) ]
    contig_results = FLYE_OUTPUT_SUMMARY.out.contig_results  // channel: [ val(meta), path(html), path(txt) ]

}
