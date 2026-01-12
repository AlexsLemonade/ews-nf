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

process metaprogram_report {
  container params.ews_nf_container
  label 'mem_24'
  time '2h'
  publishDir "${params.outdir}/metaprograms"
  input:
    path score_files
    path metaprograms_file
    path template_file
    path custom_genesets_file
    path msigdb_genesets_file
    path gene_set_category_file
  output:
    path "metaprograms-summary-report.html"
  script:
    def score_file_string = score_files.join(',')
    """
    Rscript -e "rmarkdown::render('${template_file}', \
      output_file = 'metaprograms-summary-report.html',\
      params = list(metaprograms_object_file = '${metaprograms_file}', \
                    scores_file_list = '${score_file_string}', \
                    geneset_file = '${custom_genesets_file}', \
                    msigdb_geneset_file = '${msigdb_genesets_file}', \
                    category_file = '${gene_set_category_file}'))"
    """
  stub:
    """
    touch "metaprograms-summary-report.html"
    """
}


workflow metaprograms {
  take:
    cnmf_ch  // [unique_id, sample_id, cnmf output]
    sce_ch   // [unique_id, sample_id, sce_file, annotations_file]
  main:

    // get all directories as one list
    metaprograms_input_ch = cnmf_ch.collect{ unique_id, sample_id, cnmf_output -> cnmf_output }

    // generate metaprograms
    generate_metaprograms(metaprograms_input_ch)

    // score metaprograms for each sample
    score_metaprograms(sce_ch, generate_metaprograms.out)

    // render report
    // define path to template report 
    report_template_file = "${projectDir}/templates/metaprogram-evaluation.Rmd" 
    score_ch = score_metaprograms.out.collect { unique_id, sample_id, score_files -> score_files }
    metaprogram_report(
      score_ch,
      generate_metaprograms.out,
      file(report_template_file),
      file(params.custom_genesets_file),
      file(params.msigdb_genesets_file),
      file(params.gene_set_category_file)
    )

  emit:
    metaprogram_report.out // ["metaprograms-summary-report.html"]
}
