#!/usr/bin/env nextflow

include { RUN_FLYE } from '../modules/flye'

workflow ASSEMBLE_CONTIGS {
    take:
    ch_reads_pass  // channel: [ val(meta), path(reads)]

    main:
        RUN_FLYE(ch_reads_pass)

    emit:
    flye_outputs = RUN_FLYE.out.flye_outputs  // channel: [ val(meta), path(gfa), path(fasta), path(txt) ]
}
