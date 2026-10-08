#!/usr/bin/env nextflow

include { RUN_KRAKEN2 as RUN_KRAKEN2_VIPER } from '../modules/kraken2'

workflow CLASSIFY_READS {
    take:
        ch_reads_viper   // channel: [ val(meta), path(reads), path(viper)] - reads that pass contamination check

    main:
        RUN_KRAKEN2_VIPER(ch_reads_viper,
                            'reads',
                            'viper')

    emit:
        kraken2_viper = RUN_KRAKEN2_VIPER.out.kraken2_outputs  // channel: [ val(meta), val(database_name), val(sampletype), path(kraken2 stdout) ]

}
