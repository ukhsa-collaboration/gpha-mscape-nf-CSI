#!/usr/bin/env nextflow

include { RUN_KRAKEN2 as RUN_KRAKEN2_PLUSPF } from '../modules/kraken2'
include { RUN_KRAKEN2 as RUN_KRAKEN2_VIPER } from '../modules/kraken2'

workflow CLASSIFY_CONTIGS {
    take:
    ch_contigs_pluspf // channel:[ val(meta), path(fasta), path(pluspf)]
    ch_contigs_viper // channel:[ val(meta), path(fasta), path(viper)]

    main:
        RUN_KRAKEN2_PLUSPF(ch_contigs_pluspf,
                                'contigs',
                                'pluspf')

        RUN_KRAKEN2_VIPER(ch_contigs_viper,
                            'contigs',
                            'viper')

    emit:
    kraken2_pluspf = RUN_KRAKEN2_PLUSPF.out.kraken2_outputs  // channel: [ val(meta), val(database_name), val(sampletype), path(kraken2 stdout) ]
    kraken2_viper = RUN_KRAKEN2_VIPER.out.kraken2_outputs  // channel: [ val(meta), val(database_name), val(sampletype), path(kraken2 stdout) ]

}
