#!/usr/bin/env nextflow

include { RUN_KRAKEN2 as RUN_KRAKEN2_PLUSPF_OLD } from '../modules/kraken2'
include { RUN_KRAKEN2 as RUN_KRAKEN2_PLUSPF_NEW } from '../modules/kraken2'
include { RUN_KRAKEN2 as RUN_KRAKEN2_VIPER } from '../modules/kraken2'

workflow CLASSIFY_CONTIGS {
    take:
    ch_contigs_pluspf23 // channel:[ val(meta), path(fasta), path(pluspf23)]
    ch_contigs_pluspf25 // channel:[ val(meta), path(fasta), path(pluspf25) ]
    ch_contigs_viper // channel:[ val(meta), path(fasta), path(viper)]

    main:
        RUN_KRAKEN2_PLUSPF_OLD(ch_contigs_pluspf23,
                                'contigs',
                                'pluspf2023')

        RUN_KRAKEN2_PLUSPF_NEW(ch_contigs_pluspf25,
                                'contigs',
                                'pluspf2025')

        RUN_KRAKEN2_VIPER(ch_contigs_viper,
                            'contigs',
                            'viper')

    emit:
    kraken2_pluspf23 = RUN_KRAKEN2_PLUSPF_OLD.out.kraken2_outputs  // channel: [ val(meta), val(database_name), val(sampletype), path(kraken2 stdout) ]
    kraken2_pluspf25 = RUN_KRAKEN2_PLUSPF_NEW.out.kraken2_outputs  // channel: [ val(meta), val(database_name), path(kraken2 stdout) ]
    kraken2_viper = RUN_KRAKEN2_VIPER.out.kraken2_outputs  // channel: [ val(meta), val(database_name), val(sampletype), path(kraken2 stdout) ]

}
