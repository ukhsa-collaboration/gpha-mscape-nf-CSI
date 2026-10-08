#!/usr/bin/env nextflow

include { KRAKEN2_RESULT_SUMMARY as SUMMARISE_READ_KRAKEN2_PLUSPF } from '../modules/kraken2_summary'
include { KRAKEN2_RESULT_SUMMARY as SUMMARISE_READ_KRAKEN2_VIPER } from '../modules/kraken2_summary'
include { KRAKEN2_RESULT_SUMMARY as SUMMARISE_CONTIG_KRAKEN2_PLUSPF } from '../modules/kraken2_summary'
include { KRAKEN2_RESULT_SUMMARY as SUMMARISE_CONTIG_KRAKEN2_VIPER } from '../modules/kraken2_summary'

workflow CLASSIFICATION_SUMMARY {
    take:
    ch_pluspf_reads // path to kraken2 assignments using PlusPF [ val(meta), val(database_name), val(sampletype), path(kraken2_outputs) ]
    ch_viper_reads // path to kraken2 assignments using PlusPF [ val(meta), val(database_name), val(sampletype), path(kraken2_outputs) ]

	ch_pluspf_contigs // path to kraken2 assignments using PlusPF [ val(meta), val(database_name), val(sampletype), path(kraken2_outputs) ]
    ch_viper_contigs // path to kraken2 assignments using PlusPF [ val(meta), val(database_name), val(sampletype), path(kraken2_outputs) ]

    main:
        SUMMARISE_READ_KRAKEN2_PLUSPF(ch_pluspf_reads)

        SUMMARISE_READ_KRAKEN2_VIPER(ch_viper_reads)

        SUMMARISE_CONTIG_KRAKEN2_PLUSPF(ch_pluspf_contigs)

        SUMMARISE_CONTIG_KRAKEN2_VIPER(ch_viper_contigs)

    emit:
      kraken2_read_summary_pluspf = SUMMARISE_READ_KRAKEN2_PLUSPF.out.kraken2_output_summary  // channel: [ val(meta), path(kraken2 summary) ]
      kraken2_read_summary_viper = SUMMARISE_READ_KRAKEN2_VIPER.out.kraken2_output_summary  // channel: [ val(meta), path(kraken2 summary) ]

      kraken2_contig_summary_pluspf = SUMMARISE_CONTIG_KRAKEN2_PLUSPF.out.kraken2_output_summary  // channel: [ val(meta), path(kraken2 summary) ]
      kraken2_contig_summary_viper = SUMMARISE_CONTIG_KRAKEN2_VIPER.out.kraken2_output_summary  // channel: [ val(meta), path(kraken2 summary) ]

}
