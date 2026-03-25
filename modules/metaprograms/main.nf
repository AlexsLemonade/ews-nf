#!/usr/bin/env nextflow

// Module for generating metaprograms across all samples

process generate_metaprograms {
  container params.ews_nf_container
  tag "${n_metaprograms}_metaprograms"
  label 'mem_96'
  cpus 2
  time '12h'
  publishDir "${metaprograms_publish_dir}", mode: 'copy'
  input:
    path cnmf_results_dirs
    tuple val(n_metaprograms), 
          val(metaprograms_publish_dir)
  output:
    tuple val(n_metaprograms), 
          val(metaprograms_publish_dir), 
          path("metaprograms.rds")
  script:
    def cnmf_results_string = cnmf_results_dirs.join(',')
    """
    01-generate-metaprograms.R \
      --cnmf_results_dirs ${cnmf_results_string} \
      --cnmf_k_range ${params.cnmf_k_range} \
      --n_metaprograms ${n_metaprograms} \
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
  tag "${unique_id}-${n_metaprograms}_metaprograms"
  label 'mem_24'
  cpus 2
  time '12h'
  publishDir "${metaprograms_publish_dir}/${sample_id}", mode: 'copy'
  input:
    tuple val(unique_id),
          val(sample_id),
          path(sce_file),
          path(annotations_file),
          val(n_metaprograms), 
          val(metaprograms_publish_dir), 
          path(metaprograms_file)
  output:
    tuple val(n_metaprograms),
          val(unique_id),
          val(sample_id),
          path("${unique_id}_metaprogram_scores.tsv")
  script:
    """
    02-score-metaprograms.R \
      --sce_file ${sce_file} \
      --annotations_file ${annotations_file} \
      --metaprograms_file ${metaprograms_file} \
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
  tag "${n_metaprograms}_metaprograms"
  label 'mem_24'
  time '2h'
  publishDir "${metaprograms_publish_dir}", mode: 'copy'
  input:
    tuple val(n_metaprograms), 
          val(metaprograms_publish_dir), 
          path(metaprograms_file)
    path custom_genesets_file
    path msigdb_genesets_file
    path gene_set_category_file
  output:
    tuple val(n_metaprograms), 
          val(metaprograms_publish_dir), 
          path("*metaprograms-metrics.rds")
  script:
    """
    04-metaprogram-metrics.R \
      --metaprograms_file ${metaprograms_file} \
      --sample_metadata_file ${params.sample_metadata_file} \
      --custom_geneset_file ${custom_genesets_file} \
      --msigdb_geneset_file ${msigdb_genesets_file} \
      --category_file ${gene_set_category_file} \
      --output_file "k_${n_metaprograms}_metaprograms-metrics.rds" \
      --nreps ${params.nreps} \
      --seed ${params.seed}
    """
  stub:
    """
    touch "k_${n_metaprograms}_metaprograms-metrics.rds"
    """
}

process individual_metaprogram_report {
  container params.ews_nf_container
  tag "${n_metaprograms}_metaprograms"
  label 'mem_24'
  time '2h'
  publishDir "${metaprograms_publish_dir}", mode: 'copy'
  input:
    tuple val(n_metaprograms),
          val(metaprograms_publish_dir),
          path(metaprograms_file),
          path(score_files)
    path template_file
    path custom_genesets_file
    path msigdb_genesets_file
    path gene_set_category_file
  output:
    tuple val(n_metaprograms), 
          val(metaprograms_publish_dir),
          path("metaprograms-summary-report.html")
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

process combined_metaprogram_report {
  container params.ews_nf_container
  label 'mem_24'
  time '2h'
  publishDir "${params.outdir}/metaprograms", mode: 'copy'
  input:
    path metrics_files
    path template_file
  output:
    path "combined-metaprograms-metrics.html"
  script:
    def metrics_file_string = metrics_files.join(',')
    """
    Rscript -e "rmarkdown::render('${template_file}', \
      output_file = 'combined-metaprograms-metrics.html',\
      params = list(metrics_file_list = '${metrics_file_string}', \
                    seed = '${params.seed}'))"
    """
  stub:
    """
    touch "combined-metaprograms-metrics.html"
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

    // channel of n_metaprograms and publish directory for each version of metaprograms
    n_metaprograms_ch = channel.fromList(params.n_metaprograms.toString().tokenize(","))
      .map{n_metaprograms -> 
        // define the publish directory for metaprograms
        def metaprograms_publish_dir = "${params.outdir}/metaprograms/metaprograms-${n_metaprograms}_cnmf-${cnmf_k_range_str}_${filter_str}"
        [n_metaprograms, metaprograms_publish_dir]
      }

    // get all directories as one list
    metaprograms_input_ch = cnmf_ch.collect{ unique_id, sample_id, cnmf_output -> cnmf_output }

    // generate metaprograms
    generate_metaprograms(metaprograms_input_ch, n_metaprograms_ch)

    // score metaprograms for each sample
    // make sure to score for each sample for each value of n_metaprograms
    score_input_ch = sce_ch.combine(generate_metaprograms.out)
    score_metaprograms(score_input_ch)

    // calculate metrics 
    metaprogram_metrics(
      generate_metaprograms.out, 
      file(params.custom_genesets_file),
      file(params.msigdb_genesets_file),
      file(params.gene_set_category_file)
    )

    // render report
    // define path to template reports
    individual_report_template_file = "${projectDir}/templates/metaprogram-evaluation.Rmd" 
    combined_report_template_file = "${projectDir}/templates/combined-metaprogram-metrics.Rmd"

    grouped_score_ch = score_metaprograms.out
      .map{ n_metaprograms, unique_id, sample_id, score_file -> 
        [n_metaprograms, score_file]
      }
      // group by n_metaprograms to get list of score files for each n_metaprograms
      .groupTuple(by: 0)

    report_ch = generate_metaprograms.out
      // joining on n_metaprograms
      .join(grouped_score_ch, by: 0)

    individual_metaprogram_report(
      report_ch,
      file(individual_report_template_file),
      file(params.custom_genesets_file),
      file(params.msigdb_genesets_file),
      file(params.gene_set_category_file)
    )

    // collect all metrics files for the combined report
    metrics_ch = metaprogram_metrics.out
      .collect{ n_metaprograms, metaprograms_publish_dir, metrics_file -> metrics_file }

    // render combined report
    combined_metaprogram_report(
      metrics_ch,
      file(combined_report_template_file)
    )

  emit:
    metaprograms = generate_metaprograms.out // [n_metaprograms, metaprograms_publish_dir, "metaprograms.rds"]
    metaprogram_scores = score_metaprograms.out // [n_metaprograms, unique_id, sample_id, metaprogram scores tsv]
    metaprogram_metrics = metaprogram_metrics.out // [n_metaprograms, metaprograms_publish_dir, "metaprogram-metrics.rds"]
}
