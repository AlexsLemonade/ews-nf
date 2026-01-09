#!/usr/bin/env nextflow

// Module for generating metaprograms across all samples

process generate_metaprograms {
  container params.ews_nf_container
  label 'mem_24'
  cpus 2
  time '12h'
  publishDir "${params.outdir}/metaprograms", mode: 'copy'
  input:
    path(cnmf_results_dir)
  output:
    path("metaprograms.rds")
  script:
    def cnmf_results_string = cnmf_results_dir.join(',')
    """
    01-generate-metaprograms.R \
      --cnmf_results_dir ${cnmf_results_string} \
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

process score_metaprograms {
  container params.ews_nf_container
  tag "${unique_id}"
  label 'mem_24'
  cpus 2
  time '12h'
  publishDir "${params.outdir}/metaprograms/${sample_id}", mode: 'copy'
  input:
    tuple val(unique_id),
          val(sample_id),
          path(sce_file),
          path(annotations_file)
    path metaprograms_rds
  output:
    tuple val(unique_id),
          val(sample_id),
          path("${unique_id}_metaprogram_scores.tsv")
  script:
    """
    02-score-metaprograms.R \
      --sce_file ${sce_file} \
      --annotations_file ${annotations_file} \
      --metaprograms_file ${metaprograms_rds} \
      --output_file "${unique_id}_metaprogram_scores.tsv" \
      --seed ${params.seed}
    """

  stub:
    """
    touch ${unique_id}_metaprogram_scores.tsv
    """
}


workflow metaprograms {
  take:
    cnmf_ch  // [unique_id, sample_id, cnmf output]
    sce_ch   // [unique_id, sample_id, sce_file, annotations_file]
  main:

    // get all directories as one channel
    metaprograms_input_ch = cnmf_ch.collect{ unique_id, sample_id, cnmf_output -> cnmf_output }
      .collect()

    // generate metaprograms
    generate_metaprograms(metaprograms_input_ch)

    // score metaprograms for each sample
    score_metaprograms(sce_ch, generate_metaprograms.out)

    // TODO: render report

  emit:
    generate_metaprograms.out // ["metaprograms.rds"]
}
