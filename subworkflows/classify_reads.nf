#!/usr/bin/env nextflow

include { RUN_KRAKEN2 as RUN_KRAKEN2_PLUSPF_NEW } from '../modules/kraken2'
include { RUN_KRAKEN2 as RUN_KRAKEN2_VIPER } from '../modules/kraken2'

workflow CLASSIFY_READS {
    take:
        ch_reads_pass   // channel: [ val(meta), path(reads)] - reads that pass contamination check
        ch_pluspf25 // path to latest PlusPF database
        ch_viper // path to curated viral database

    main:
        RUN_KRAKEN2_PLUSPF_NEW(ch_reads_pass,
                                ch_pluspf25,
                                'reads',
                                'pluspf2025')

        RUN_KRAKEN2_VIPER(ch_reads_pass,
                            ch_viper,
                            'reads',
                            'viper')

    emit:
        kraken2_pluspf25 = RUN_KRAKEN2_PLUSPF_NEW.out.kraken2_outputs  // channel: [ val(meta), val(database_name), val(sampletype), path(kraken2 stdout) ]
        kraken2_viper = RUN_KRAKEN2_VIPER.out.kraken2_outputs  // channel: [ val(meta), val(database_name), val(sampletype), path(kraken2 stdout) ]

}
