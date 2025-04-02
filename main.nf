#!/usr/bin/env nextflow

// **** Parameter checks ****
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
param_error = false

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

if (!params.scpca_nf_dir.exists()) {
  log.error "scpca-nf output directory does not exist: ${scpca_nf_dir}"
  param_error = true
}

if (param_error) {
  System.exit(1)
}

// stub test
process test_setup {
  publishDir "${params.outdir}/${sample_id}"
  input:
    tuple val(sample_id), path(library_files)
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

  // create a channel of [sample id, [path to processed SCE files]]
  sample_ch = Channel.fromList(Utils.getLibraryFiles(params.scpca_nf_dir))
    .filter{
      run_all || (it[0] in sample_ids)
    }

  // TODO: Remove this when we add our first module
  // run stub test to make sure setup is good
  test_setup(sample_ch)

}
