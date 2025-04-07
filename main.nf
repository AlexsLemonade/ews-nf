#!/usr/bin/env nextflow

// **** Parameter checks ****
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
param_error = false

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

// TODO: remove this once we add in a real process
// for now just use this to test the configuration 
process test_setup {
  publishDir "${params.outdir}/${sample_id}"
  input:
    tuple val(sample_id), val(library_id), path(library_file)
  output:
    path test_output
  script: 
    test_output = "test.txt"
    """
    touch ${test_output}
    """
  stub:
    test_output = "test.txt"
    """
    touch ${test_output}
    """
}


workflow {
  sample_ids = params.sample_ids?.tokenize(';, ') ?: []
  run_all = sample_ids[0] == "All"

  println "scpca_nf_dir = ${params.scpca_nf_dir}"

  // create a channel of [sample id, library_id, path to processed SCE file]
  library_ch = Channel.fromPath("${params.scpca_nf_dir}/**/**_processed.rds", type: 'file')
    .map{ file -> tuple(
      file.parent.name, // extract sample id from file path
      file.name.replace("_processed.rds", ""), // extract library id
      file // full path to file
    )}

  // TODO: Remove this when we add our first module
  // run stub test to make sure setup is good
  test_setup(library_ch)

}
