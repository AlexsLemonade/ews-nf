#!/usr/bin/env nextflow

// **** Included workflows and processes from modules ****
include { aucell } from './modules/aucell'
include { build_ref } from './modules/build-reference'
include { infercnv } from './modules/infercnv'
include { summary_report } from './modules/summary-report'

// **** Parameter checks ****
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
param_error = false

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)


workflow {
  sample_ids = params.sample_ids?.tokenize(';, ') ?: []
  run_all = sample_ids[0] == "All"

  // create a channel of [unique_id, sample id, path to processed SCE file]
  library_ch = Channel.fromPath("${params.scpca_nf_dir}/**/**_processed.rds", type: 'file')
    .map{ file -> tuple(
      file.name.replace("_processed.rds", ""), // extract unique id which will be either library_id or <library_id>-<sample_id>
      file.parent.name, // extract sample id from file path
      file // full path to file
    )}
    .filter{ run_all || (it[0] in sample_ids) } // filter to only specified samples

  // build reference of normal cells
  build_ref(library_ch)

  // run aucell 
  aucell(library_ch, build_ref.out)
  
  // run infercnv 
  infercnv(library_ch, build_ref.out)

  // generate report 
  summary_report(library_ch, aucell.out.query, aucell.out.ref)
  
}
