#!/usr/bin/env nextflow

include { WRITE_CONFIG } from '../modules/report_config'

workflow WRITE_REPORT {
    take:
    ch_config_input 

    main:
        WRITE_CONFIG(ch_config_input)


}
