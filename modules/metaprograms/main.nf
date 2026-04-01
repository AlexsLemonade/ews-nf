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
          val(metaprograms_publish_dir),
          val(unique_id),
          val(sample_id),
          path("${unique_id}_metaprogram_scores.tsv"),
          path("${unique_id}_scores_background.tsv")
  script:
    """
    02-score-metaprograms.R \
      --sce_file ${sce_file} \
      --annotations_file ${annotations_file} \
      --metaprograms_file ${metaprograms_file} \
      --output_file "${unique_id}_metaprogram_scores.tsv" \
      --seed ${params.seed}

    03-calculate-score-background.R \
      --sce_file ${sce_file} \
      --annotations_file ${annotations_file} \
      --metaprograms_file ${metaprograms_file} \
      --output_file "${unique_id}_scores_background.tsv" \
      --seed ${params.seed}
    """

  stub:
    """
    touch ${unique_id}_metaprogram_scores.tsv
    touch ${unique_id}_scores_background.tsv
    """
}

process combine_score_files {
  container params.ews_nf_container
  tag "${n_metaprograms}_metaprograms"
  publishDir "${metaprograms_publish_dir}", mode: 'copy'
  label 'mem_24'
  cpus 2
  time '12h'
  input:
    tuple val(n_metaprograms), 
          val(metaprograms_publish_dir),
          path(score_files), // list of score files for each sample, 
          path(background_score_files) // list of background score files for each sample
  output:
    tuple val(n_metaprograms),
          path("k-${n_metaprograms}_metaprogram_scores.tsv"),
          path("k-${n_metaprograms}_background_scores.tsv")
  script:
    """
    03a-combine-score-files.R \
      --input_tsv_files ${score_files.join(',')} \
      --output_file "k-${n_metaprograms}_metaprogram_scores.tsv"

    03a-combine-score-files.R \
      --input_tsv_files ${background_score_files.join(',')} \
      --output_file "k-${n_metaprograms}_background_scores.tsv"
    """

  stub:
    """
    touch "k-${n_metaprograms}_metaprogram_scores.tsv"
    touch "k-${n_metaprograms}_background_scores.tsv"
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
          path("k-${n_metaprograms}_metaprograms-metrics.rds")
  script:
    """
    04-metaprogram-metrics.R \
      --metaprograms_file ${metaprograms_file} \
      --sample_metadata_file ${params.sample_metadata_file} \
      --custom_geneset_file ${custom_genesets_file} \
      --msigdb_geneset_file ${msigdb_genesets_file} \
      --category_file ${gene_set_category_file} \
      --output_file "k-${n_metaprograms}_metaprograms-metrics.rds" \
      --nreps ${params.nreps} \
      --seed ${params.seed}
    """
  stub:
    """
    touch "k-${n_metaprograms}_metaprograms-metrics.rds"
    """
}

// used for comparing across k values
process combined_metaprogram_report {
  container params.ews_nf_container
  label 'mem_24'
  time '2h'
  publishDir "${params.outdir}/metaprograms", mode: 'copy'
  input:
    path metrics_files
    path scores_files
    path background_scores_files
    path template_file
  output:
    path "combined-metaprograms-metrics.html", emit: combined_report
    env "OPTIMAL_K_VALUE", emit: optimal_k
  script:
    def metrics_file_string = metrics_files.join(',')
    def scores_file_string = scores_files.join(',')
    def background_scores_files_string = background_scores_files.join(',')
    """
    Rscript -e "rmarkdown::render('${template_file}', \
      output_file = 'combined-metaprograms-metrics.html',\
      params = list(metrics_file_list = '${metrics_file_string}', \
                    scores_file_list = '${scores_file_string}', \
                    background_scores_file_list = '${background_scores_files_string}', \
                    sample_metadata_file = '${params.sample_metadata_file}', \
                    seed = '${params.seed}'))"

    # save optimal k as a variable to be used to score cells 
    OPTIMAL_K_VALUE=\$(cat "optimal_k_value.txt")
    """
  stub:
    """
    touch "combined-metaprograms-metrics.html"
    OPTIMAL_K_VALUE=10
    """
}

// used for individual k values
process individual_metaprogram_report {
  container params.ews_nf_container
  tag "${n_metaprograms}_metaprograms"
  label 'mem_24'
  time '2h'
  publishDir "${params.outdir}/metaprograms", mode: 'copy'
  input:
    tuple val(n_metaprograms),
          path(metaprograms_file),
          path(score_files)
    path template_file
    path custom_genesets_file
    path msigdb_genesets_file
    path gene_set_category_file
  output:
    tuple val(n_metaprograms), 
          path("k-${n_metaprograms}_metaprograms-summary-report.html")
  script:
    def score_file_string = score_files.join(',')
    """
    Rscript -e "rmarkdown::render('${template_file}', \
      output_file = 'k-${n_metaprograms}_metaprograms-summary-report.html',\
      params = list(metaprograms_object_file = '${metaprograms_file}', \
                    scores_file_list = '${score_file_string}', \
                    geneset_file = '${custom_genesets_file}', \
                    msigdb_geneset_file = '${msigdb_genesets_file}', \
                    category_file = '${gene_set_category_file}'))"
    """
  stub:
    """
    touch "k-${n_metaprograms}_metaprograms-summary-report.html"
    """
}

// score and annotate objects 
process score_and_annotate_sce {
  container params.ews_nf_container
  tag "${unique_id}"
  label 'mem_24'
  time '2h'
  publishDir "${params.outdir}/metaprograms/${sample_id}", mode: 'copy'
  input:
    tuple val(unique_id),
          val(sample_id),
          path(sce_file),
          path(annotations_file),
          path(score_file)
  output:
    path("${unique_id}_annotated.sce")
  script:
    """
    05-score-and-annotate-sce.R \
      --sce_file ${sce_file} \
      --annotations_file ${annotations_file} \
      --score_file ${score_file} \
      --output_file "${unique_id}_annotated.sce"
    """
  stub:
    """
    touch "${unique_id}_annotated.sce"
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

    grouped_score_ch = score_metaprograms.out
      .map{ n_metaprograms, metaprograms_publish_dir, unique_id, sample_id, score_file, background_score_file -> 
        [n_metaprograms, metaprograms_publish_dir, score_file, background_score_file]
      }
      // group by n_metaprograms to get list of score files and background score files for each n_metaprograms
      .groupTuple(by: [0,1]) // group by n_metaprograms and publish dir


    // create a single data frame for all observed and all background scores
    combine_score_files(grouped_score_ch)

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

    // collect all metrics files for the combined report
    metrics_ch = metaprogram_metrics.out
      .collect{ n_metaprograms, metaprograms_publish_dir, metrics_file -> metrics_file }

    // collect all score and background score files for the report
    scores_ch = combine_score_files.out
      .collect{ n_metaprograms, score_file, background_score_file -> score_file }

    background_scores_ch = combine_score_files.out
      .collect{ n_metaprograms, score_file, background_score_file -> background_score_file }

    // render combined report
    combined_metaprogram_report(
      metrics_ch,
      scores_ch, 
      background_scores_ch,
      file(combined_report_template_file)
    )

    // make a channel of sce file, annotation file and score file 
    // only for the optimal k value 
    score_and_annotate_ch = grouped_score_ch
      // drop publish directory and background score files 
      .map{ n_metaprograms, metaprograms_publish_dir, score_files, background_score_files -> 
        [n_metaprograms, score_files]
      }
      .combine(combined_metaprogram_report.out.optimal_k) // combine with optimal k value from combined report
      // filter to only get score outputs for optimal k value
      .filter{ n_metaprograms, score_files, optimal_k -> n_metaprograms == optimal_k }
      .map{ it.dropRight(1) } // drop extra k 
      // transpose channel to split up score files [[n_metaprograms, score_file], [n_metaprograms, score_file]]
      .transpose()
      // pull out unique id from score file name and use that to join with sce_ch 
      .map{ n_metaprograms, score_file -> 
        def unique_id = score_file.name.split("_")[0] // extract unique id from score file name
        [unique_id, score_file]
      }
      // join with sce channel to add in sce file and annotation file for each sample
      .join(sce_ch, by: 0)
      .map{ unique_id, score_file, sample_id, sce_file, annotations_file -> 
        [unique_id, sample_id, sce_file, annotations_file, score_file]
      }

    // output sces with annotations and scores for the optimal k value
    score_and_annotate_sce(score_and_annotate_ch)

    // only render a report for the optimal k value
    report_ch = generate_metaprograms.out
      .combine(combined_metaprogram_report.out.optimal_k)
      .filter { n_metaprograms, metaprograms_publish_dir, metaprograms_file, optimal_k -> n_metaprograms == optimal_k }
      // joining on n_metaprograms
      // drop the background files before joining
      .join(grouped_score_ch.map { it.dropRight(1) }, by: [0, 1])
      // drop publish directory and extra value of k 
      // instead we ouptut to the main metaprograms directory
      .map{ n_metaprograms, metaprograms_publish_dir, metaprograms_file, optimal_k, score_files -> 
        [n_metaprograms, metaprograms_file, score_files]
      }

    individual_metaprogram_report(
      report_ch,
      file(individual_report_template_file),
      file(params.custom_genesets_file),
      file(params.msigdb_genesets_file),
      file(params.gene_set_category_file)
    )

  emit:
    metaprograms = generate_metaprograms.out // [n_metaprograms, metaprograms_publish_dir, "metaprograms.rds"]
    metaprogram_scores = score_metaprograms.out // [n_metaprograms, unique_id, sample_id, metaprogram scores tsv, background scores tsv]
    metaprogram_metrics = metaprogram_metrics.out // [n_metaprograms, metaprograms_publish_dir, "metaprogram-metrics.rds"]
}
