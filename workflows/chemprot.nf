#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CHEMPROT_TPP } from "../modules/local/chemprot_tpp/main"
include { CHEMPROT_DOSE } from "../modules/local/chemprot_dose/main"
include { CHEMPROT_ABPP } from "../modules/local/chemprot_abpp/main"

workflow CHEMPROT {
    take:
    ch_msstats

    main:
    if (params.analysis_type == "tpp") {
        CHEMPROT_TPP(ch_msstats, params.tpp_temperatures, params.tpp_min_r2, params.min_peptides)
    } else if (params.analysis_type == "dose_response") {
        CHEMPROT_DOSE(ch_msstats, params.compound_conc, params.dose_response_min_r2, params.min_peptides)
    } else if (params.analysis_type == "abpp") {
        CHEMPROT_ABPP(ch_msstats, params.abpp_condition_col, params.fdr_threshold)
    }
}
