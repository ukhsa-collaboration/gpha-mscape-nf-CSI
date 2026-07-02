include { UNCLASSIFIED_INVESTIGATION } from './workflows/unclassified'

workflow {
    if (!params.input) {
        error "Missing required parameter: --input"
    }

    if (!params.pluspf) {
        error "Missing required parameter: --pluspf"
    }

    if (!params.pluspf_recent) {
        error "Missing required parameter: --pluspf_recent"
    }

    if (!params.viper) {
        error "Missing required parameter: --viper"
    }


    def ch_samples = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> tuple([id: row.climb_id], row) }

    def ch_sample_inputs = ch_samples.multiMap { meta, row ->
        reads: [meta, file(row.fastq)]
        scylla_output: [meta, 'pluspf2023', 'reads', file(row.kraken_stdout)]
    }

    def ch_ref_fasta = Channel.fromPath(params.ref_fasta)

    def ch_pluspf23 = Channel.fromPath(params.pluspf)

    def ch_pluspf25 = Channel.fromPath(params.pluspf_recent)

    def ch_viper = Channel.fromPath(params.viper)

    UNCLASSIFIED_INVESTIGATION(
        ch_sample_inputs.reads,
        ch_ref_fasta,
        ch_pluspf23,
        ch_pluspf25,
        ch_viper,
        ch_sample_inputs.scylla_output
    )
}
