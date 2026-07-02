#!/usr/bin/env nextflow

include { SANKEY } from '../modules/sankey'

workflow PLOT_CLASSIFICATION_RESULTS {
    take:
    ch_mapping_stats
    ch_read2contig

    ch_pluspf23_reads // path to kraken2 assignments using PlusPF [ val(meta), val(database_name), val(sampletype), path(kraken2_outputs) ]
    ch_pluspf25_reads // path to kraken2 assignments using latest PlusPF database  [ val(meta), val(sampletype), val(database_name), path(kraken2_outputs) ]
    ch_viper_reads // path to kraken2 assignments using PlusPF [ val(meta), val(database_name), val(sampletype), path(kraken2_outputs) ]

	ch_pluspf23_contigs // path to kraken2 assignments using PlusPF [ val(meta), val(database_name), val(sampletype), path(kraken2_outputs) ]
    ch_pluspf25_contigs // path to kraken2 assignments using latest PlusPF database  [ val(meta), val(sampletype), val(database_name), path(kraken2_outputs) ]
    ch_viper_contigs // path to kraken2 assignments using PlusPF [ val(meta), val(database_name), val(sampletype), path(kraken2_outputs) ]

    main:
        SANKEY(ch_mapping_stats,
               ch_read2contig,
               ch_pluspf23_reads,
               ch_pluspf25_reads,
               ch_viper_reads,
               ch_pluspf23_contigs,
               ch_pluspf25_contigs,
               ch_viper_contigs)

}
