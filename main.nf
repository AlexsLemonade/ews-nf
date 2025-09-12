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
    .filter{ run_all || (it[1] in sample_ids) } // filter to only specified samples
    .branch{
      run_workflow: it[2].size() > 0 || it[2].name.startsWith("STUB")
      skip_workflow: true
    }

  // send unique ids for all libraries that don't have enough cells to log 
  library_ch.skip_workflow
    .subscribe{
      log.error("There are no cells found in the processed object for ${it[0]}.")
    }

  run_ch = library_ch.run_workflow

  // build reference of normal cells
  build_ref(run_ch)

  // run aucell 
  aucell(run_ch, build_ref.out)
  
  // run infercnv 
  infercnv(run_ch, build_ref.out)

  // generate report 
  summary_report(run_ch, aucell.out.query, aucell.out.ref, infercnv.out, build_ref.out)
  
}
