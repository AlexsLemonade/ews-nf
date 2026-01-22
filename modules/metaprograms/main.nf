#!/usr/bin/env nextflow

// Module for generating metaprograms across all samples

process generate_metaprograms {
  container params.ews_nf_container
  label 'mem_96'
  cpus 2
  time '12h'
  publishDir "${metaprograms_publish_dir}", mode: 'copy'
  input:
    path cnmf_results_dirs
    val metaprograms_publish_dir
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
  publishDir "${metaprograms_publish_dir}/${sample_id}", mode: 'copy'
  input:
    tuple val(unique_id),
          val(sample_id),
          path(sce_file),
          path(annotations_file)
    path metaprograms_rds
    val metaprograms_publish_dir
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

process metaprogram_metrics {
  container params.ews_nf_container
  label 'mem_24'
  time '2h'
  publishDir "${metaprograms_publish_dir}", mode: 'copy'
  input:
    path metaprograms_file
    path custom_genesets_file
    path msigdb_genesets_file
    path gene_set_category_file
    val metaprograms_publish_dir
  output:
    path "metaprograms-metrics.rds"
  script:
    """
    03-metaprogram-metrics.R \
      --metaprograms_file ${metaprograms_file} \
      --custom_geneset_file ${custom_genesets_file} \
      --msigdb_geneset_file ${msigdb_genesets_file} \
      --category_file ${gene_set_category_file} \
      --output_file "metaprograms-metrics.rds" \
      --seed ${params.seed}
    """
  stub:
    """
    touch "metaprograms-metrics.rds"
    """
}

process metaprogram_report {
  container params.ews_nf_container
  label 'mem_24'
  time '2h'
  publishDir "${metaprograms_publish_dir}", mode: 'copy'
  input:
    path score_files
    path metaprograms_file
    path template_file
    path custom_genesets_file
    path msigdb_genesets_file
    path gene_set_category_file
    val metaprograms_publish_dir
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

    // prepare strings for nanming the metaprogram output directory
    cnmf_k_range_str = params.cnmf_k_range.replaceAll(",", "-")
    filter_str = params.filter_spectra ? "filtered" : "unfiltered"

    // define the publish directory for metaprograms
    metaprograms_publish_dir = "${params.outdir}/metaprograms/metaprograms-${params.n_metaprograms}_cnmf-${cnmf_k_range_str}_${filter_str}"

    // get all directories as one list
    metaprograms_input_ch = cnmf_ch.collect{ unique_id, sample_id, cnmf_output -> cnmf_output }

    // generate metaprograms
    generate_metaprograms(metaprograms_input_ch, metaprograms_publish_dir)

    // score metaprograms for each sample
    score_metaprograms(sce_ch, generate_metaprograms.out, metaprograms_publish_dir)

    // calculate metrics 
    metaprogram_metrics(
      generate_metaprograms.out, 
      file(params.custom_genesets_file),
      file(params.msigdb_genesets_file),
      file(params.gene_set_category_file),
      metaprograms_publish_dir
    )

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
      file(params.gene_set_category_file),
      metaprograms_publish_dir
    )

  emit:
    metaprograms = generate_metaprograms.out // ["metaprograms.rds"]
    metaprogram_scores = score_metaprograms.out // [unique_id, sample_id, metaprogram scores tsv]
    metaprogram_metrics = metaprogram_metrics.out // ["metaprogram-metrics.rds"]
}
