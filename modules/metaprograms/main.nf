#!/usr/bin/env nextflow

// Module for generating metaprograms across all samples

process generate_metaprograms {
  container params.ews_nf_container
  label 'mem_96'
  cpus 2
  time '12h'
  publishDir "${params.outdir}/metaprograms", mode: 'copy'
  input:
    path cnmf_results_dirs
  output:
    path "metaprograms.rds"
  script:
    def cnmf_results_string = cnmf_results_dirs.join(',')
    """
    01-generate-metaprograms.R \
      --cnmf_results_dirs ${cnmf_results_string} \
      --cnmf_k_range ${params.cnmf_k_range} \
      --n_metaprograms ${params.n_metaprograms} \
      --filter_spectra ${params.filter_spectra} \
      --orphan_cutoff ${params.orphan_cutoff} \
      --output_file "metaprograms.rds" \
      --seed ${params.seed}
    """

  stub:
    """
    touch metaprograms.rds
    """
}


workflow metaprograms {
  take:
    cnmf_ch  // [unique_id, sample_id, cnmf output]
  main:

    // get all directories as one list
    metaprograms_input_ch = cnmf_ch.collect{ unique_id, sample_id, cnmf_output -> cnmf_output }

    // generate metaprograms
    generate_metaprograms(metaprograms_input_ch)

    // TODO: render report

  emit:
    generate_metaprograms.out // ["metaprograms.rds"]
}
