#!/usr/bin/env nextflow

include { READ_MAPPING } from '../subworkflows/read_mapping'
include { TOBACCO_CONTAMINATION_CHECK } from '../modules/tobacco_check'
include { CLASSIFY_READS } from '../subworkflows/classify_reads'
include { ASSEMBLE_CONTIGS } from '../subworkflows/assemble_contigs'
include { CONTIG_SUMMARY } from '../subworkflows/contig_summary'
include { CLASSIFY_CONTIGS } from '../subworkflows/classify_contigs'
include { CLASSIFICATION_SUMMARY } from '../subworkflows/classification_summary'
include { PLOT_CLASSIFICATION_RESULTS } from '../subworkflows/plot_classification_results'

workflow UNCLASSIFIED_INVESTIGATION {
    take:
    ch_reads   // channel: [ val(meta), path(reads) ]
    ch_ref_fasta  // reference fasta path
    ch_pluspf23 // path to old PlusPF database
    ch_pluspf25 // path to latest PlusPF database
    ch_viper // path to curated viral database
    ch_classified_reads_pluspf23 // path to scylla kraken2 assignments

    main:
        READ_MAPPING(ch_reads,
                     ch_ref_fasta)

        TOBACCO_CONTAMINATION_CHECK(READ_MAPPING.out.mappy_results,
                                    ch_reads)

        TOBACCO_CONTAMINATION_CHECK.out.tobacco_contamination_results
            .branch { meta, reads, tobacco ->
               cont: tobacco.text.contains("False")
                   return tuple ( meta, reads)
               end: tobacco.text.contains("True")
                   return tuple ( meta, reads)
            }
            .set{ ch_pipeline_status }

        // will only run if unclassified reads are not tobacco contamination
        CLASSIFY_READS(ch_pipeline_status.cont,
                       ch_pluspf25,
                       ch_viper)

        ASSEMBLE_CONTIGS(ch_pipeline_status.cont)

        CONTIG_SUMMARY(ch_pipeline_status.cont,
                       ASSEMBLE_CONTIGS.out.flye_outputs)

        ch_contigs = ASSEMBLE_CONTIGS.out.flye_outputs
            .map { meta, gfa, assembly, info ->
                tuple ( meta, assembly )
            }

        CLASSIFY_CONTIGS(ch_contigs,
                         ch_pluspf23,
                         ch_pluspf25,
                         ch_viper)

        CLASSIFICATION_SUMMARY(ch_classified_reads_pluspf23,
                               CLASSIFY_READS.out.kraken2_pluspf25,
                               CLASSIFY_READS.out.kraken2_viper,
                               CLASSIFY_CONTIGS.out.kraken2_pluspf23,
                               CLASSIFY_CONTIGS.out.kraken2_pluspf25,
                               CLASSIFY_CONTIGS.out.kraken2_viper)

        PLOT_CLASSIFICATION_RESULTS(READ_MAPPING.out.mappy_stats,
                                    CONTIG_SUMMARY.out.read2contigs,
                                    CLASSIFICATION_SUMMARY.out.kraken2_read_summary_pluspf23,
                                    CLASSIFICATION_SUMMARY.out.kraken2_read_summary_pluspf25,
                                    CLASSIFICATION_SUMMARY.out.kraken2_read_summary_viper,
                                    CLASSIFICATION_SUMMARY.out.kraken2_contig_summary_pluspf23,
                                    CLASSIFICATION_SUMMARY.out.kraken2_contig_summary_pluspf25,
                                    CLASSIFICATION_SUMMARY.out.kraken2_contig_summary_viper)

}
