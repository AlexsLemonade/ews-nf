#!/usr/bin/env nextflow

// Module for comparing metaprograms generated using all samples vs. subsets of samples

process comparison_metrics {
  container params.ews_nf_container
  publishDir "${params.outdir}/metaprograms/${subset_id}/", mode: 'copy', overwrite: true
  tag "${subset_id}"
  label 'mem_16'
  cpus 2
  time '12h'
  input:
    tuple val(subset_id), 
          path(subset_metaprograms_file, stageAs: "subset_file"),
          path(full_metaprograms_file, stageAs: "full_file")
  output:
    tuple val(subset_id),
          path("${subset_id}_metaprograms-comparison-metrics.tsv")
  script:
    """
    01-compare-metaprograms.R \
      --subset_id ${subset_id} \
      --subset_metaprograms_file ${subset_metaprograms_file} \
      --full_metaprograms_file ${full_metaprograms_file} \
      --output_file "${subset_id}_metaprograms-comparison-metrics.tsv"
    """

  stub:
    """
    touch "${subset_id}_metaprograms-comparison-metrics.tsv"
    """
}

workflow compare_metaprograms {
  take:
    subset_ch  // [subset id, subset metaprograms file, full metaprograms file]
  main:

    // generate a file comparing each subset mp to all full mps 
    comparison_metrics(subset_ch)

  emit:
    comparison_metrics.out

}