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
        ch_mapping_input = ch_reads
            .combine(ch_ref_fasta)

        READ_MAPPING(ch_mapping_input)

        TOBACCO_CONTAMINATION_CHECK(READ_MAPPING.out.mappy_results,
                                    ch_reads)


        TOBACCO_CONTAMINATION_CHECK.out.tobacco_contamination_results
            .branch { meta, reads, tobacco ->
               cont: tobacco.text.contains("False")
                   return tuple ( meta, reads )
               end: tobacco.text.contains("True")
                   return tuple ( meta, reads )
            }
            .set{ ch_pipeline_status }

        ch_pipeline_status.cont.view { it -> "cont value: ${it}" }
        ch_pipeline_status.end.view { it -> "end value: ${it}" }

        // make a channel for samples that pass with database otherwise does not run per sample
        ch_read_pass_pluspf25 = ch_pipeline_status.cont
            .combine(ch_pluspf25)

        // make a channel for samples that pass with database
        ch_read_pass_viper = ch_pipeline_status.cont
            .combine(ch_viper)


        // will only run if unclassified reads are not tobacco contamination
        CLASSIFY_READS(ch_read_pass_pluspf25,
                       ch_read_pass_viper)

        ASSEMBLE_CONTIGS(ch_pipeline_status.cont)

        CONTIG_SUMMARY(ch_pipeline_status.cont,
                       ASSEMBLE_CONTIGS.out.flye_outputs)

       ch_contigs = ASSEMBLE_CONTIGS.out.flye_outputs
           .map { meta, gfa, assembly, info ->
               tuple ( meta, assembly )
           }

       ch_contigs_pluspf23 = ch_contigs
           .combine(ch_pluspf23)

       ch_contigs_pluspf25 = ch_contigs
           .combine(ch_pluspf25)

       ch_contigs_viper = ch_contigs
           .combine(ch_viper)

       CLASSIFY_CONTIGS(ch_contigs_pluspf23,
                        ch_contigs_pluspf25,
                        ch_contigs_viper)


      CLASSIFICATION_SUMMARY(ch_classified_reads_pluspf23,
                             CLASSIFY_READS.out.kraken2_pluspf25,
                             CLASSIFY_READS.out.kraken2_viper,
                             CLASSIFY_CONTIGS.out.kraken2_pluspf23,
                             CLASSIFY_CONTIGS.out.kraken2_pluspf25,
                             CLASSIFY_CONTIGS.out.kraken2_viper)

      ch_sankey_input = READ_MAPPING.out.mappy_stats
          .join(CONTIG_SUMMARY.out.read2contigs)
          .join(CLASSIFICATION_SUMMARY.out.kraken2_read_summary_pluspf23)
          .join(CLASSIFICATION_SUMMARY.out.kraken2_read_summary_pluspf25)
          .join(CLASSIFICATION_SUMMARY.out.kraken2_read_summary_viper)
          .join(CLASSIFICATION_SUMMARY.out.kraken2_contig_summary_pluspf23)
          .join(CLASSIFICATION_SUMMARY.out.kraken2_contig_summary_pluspf25)
          .join(CLASSIFICATION_SUMMARY.out.kraken2_contig_summary_viper)
          .view()


      PLOT_CLASSIFICATION_RESULTS(ch_sankey_input)

}
