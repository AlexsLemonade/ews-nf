#!/usr/bin/env nextflow

// Module for generating metaprograms across all samples

process generate_metaprograms {
  container params.ews_nf_container
  tag "${subset_id}-${n_metaprograms}_metaprograms"
  label 'mem_96'
  cpus 2
  time '12h'
  publishDir "${metaprograms_publish_dir}", mode: 'copy', enabled: !params.subset_mode
  input:
    tuple val(subset_id), 
          val(n_metaprograms), 
          val(unique_id_tuple),
          path(cnmf_output_dirs),
          val(metaprograms_publish_dir)
  output:
    tuple val(subset_id),
          val(n_metaprograms), 
          val(unique_id_tuple),
          val(metaprograms_publish_dir), 
          path("k-${n_metaprograms}_metaprograms.rds")
  script:
    def cnmf_results_string = cnmf_output_dirs.join(',')
    """
    01-generate-metaprograms.R \
      --cnmf_results_dirs ${cnmf_results_string} \
      --cnmf_k_range ${params.cnmf_k_range} \
      --n_metaprograms ${n_metaprograms} \
      --filter_spectra ${params.filter_spectra} \
      --orphan_cutoff ${params.orphan_cutoff} \
      --output_file "k-${n_metaprograms}_metaprograms.rds" \
      --seed ${params.seed}
    """

  stub:
    """
    touch "k-${n_metaprograms}_metaprograms.rds"
    """
}

process score_metaprograms {
  container params.ews_nf_container
  tag "${subset_id}-${unique_id}-${n_metaprograms}_metaprograms"
  label 'mem_24'
  cpus 2
  time '12h'
  input:
    tuple val(unique_id), 
          val(subset_id), 
          val(n_metaprograms),  
          val(metaprograms_publish_dir),
          path(metaprograms_file),
          val(sample_id),
          path(sce_file), 
          path(annotations_file)
  output:
    tuple val(subset_id), 
          val(n_metaprograms),
          val(metaprograms_publish_dir),
          val(unique_id),
          val(sample_id),
          path("${unique_id}_metaprogram_scores.tsv.gz"),
          path("${unique_id}_scores_background.tsv.gz")
  script:
    """
    02-score-metaprograms.R \
      --sce_file ${sce_file} \
      --annotations_file ${annotations_file} \
      --unique_id ${unique_id} \
      --metaprograms_file ${metaprograms_file} \
      --output_file "${unique_id}_metaprogram_scores.tsv.gz" \
      --seed ${params.seed}

    03-calculate-score-background.R \
      --sce_file ${sce_file} \
      --annotations_file ${annotations_file} \
      --unique_id ${unique_id} \
      --metaprograms_file ${metaprograms_file} \
      --output_file "${unique_id}_scores_background.tsv.gz" \
      --seed ${params.seed}
    """

  stub:
    """
    touch ${unique_id}_metaprogram_scores.tsv.gz
    touch ${unique_id}_scores_background.tsv.gz
    """
}

process combine_score_files {
  container params.ews_nf_container
  tag "${subset_id}-${n_metaprograms}_metaprograms"
  publishDir "${metaprograms_publish_dir}", mode: 'copy', enabled: !params.subset_mode
  label 'mem_24'
  cpus 2
  time '12h'
  input:
    tuple val(subset_id),
          val(n_metaprograms), 
          val(metaprograms_publish_dir),
          path(score_files), // list of score files for each sample, 
          path(background_score_files) // list of background score files for each sample
  output:
    tuple val(subset_id),
          val(n_metaprograms),
          path("k-${n_metaprograms}_metaprogram_scores.tsv.gz"),
          path("k-${n_metaprograms}_background_scores.tsv.gz")
  script:
    """
    03a-combine-score-files.R \
      --input_tsv_files ${score_files.join(',')} \
      --output_file "k-${n_metaprograms}_metaprogram_scores.tsv.gz"

    03a-combine-score-files.R \
      --input_tsv_files ${background_score_files.join(',')} \
      --output_file "k-${n_metaprograms}_background_scores.tsv.gz"
    """

  stub:
    """
    touch "k-${n_metaprograms}_metaprogram_scores.tsv.gz"
    touch "k-${n_metaprograms}_background_scores.tsv.gz"
    """
}



process metaprogram_metrics {
  container params.ews_nf_container
  tag "${subset_id}-${n_metaprograms}_metaprograms"
  label 'mem_24'
  time '2h'
  publishDir "${metaprograms_publish_dir}", mode: 'copy', enabled: !params.subset_mode
  input:
    tuple val(subset_id),
          val(n_metaprograms), 
          val(unique_id_tuple),
          val(metaprograms_publish_dir), 
          path(metaprograms_file)
    path custom_genesets_file
    path msigdb_genesets_file
    path gene_set_category_file
  output:
    tuple val(subset_id),
          val(n_metaprograms),
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
  tag "${subset_id}"
  publishDir "${params.outdir}/metaprograms/${subset_id}", mode: 'copy'
  input:
    tuple val(subset_id),
          val(n_metaprograms),
          path(metrics_files), // list of metrics files for each n_metaprograms
          path(scores_files), // list of combined score files for each n_metaprograms
          path(background_scores_files) // list of combined background score files for each n_metaprograms
    path template_file
  output:  
    tuple val(subset_id),
          env("OPTIMAL_K_VALUE"), 
          path("metrics-heatmap.png")
          path("combined-metaprograms-metrics.html")
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
                    output_heatmap_file = 'metrics-heatmap.png', \
                    output_k_file = 'optimal_k_value.txt', \
                    seed = '${params.seed}'))"

    # save optimal k as a variable to be used to score cells 
    OPTIMAL_K_VALUE=\$(cat "optimal_k_value.txt")
    """
  stub:
    """
    touch "combined-metaprograms-metrics.html"
    touch "metrics-heatmap.png"
    OPTIMAL_K_VALUE=10
    """
}

// used for individual k values
process individual_metaprogram_report {
  container params.ews_nf_container
  tag "${subset_id}-${n_metaprograms}_metaprograms"
  label 'mem_24'
  time '2h'
  publishDir "${params.outdir}/metaprograms/${subset_id}", mode: 'copy'
  input:
    tuple val(subset_id), 
          val(n_metaprograms),
          path(metaprograms_file),
          path(metrics_file),
          path(combined_score_file)
    path template_file
    path gene_set_category_file
  output:
    tuple val(subset_id),
          val(n_metaprograms), 
          path(metaprograms_file),
          path(metrics_file),
          path(combined_score_file),
          path("k-${n_metaprograms}_metaprograms-summary-report.html")
  script:
    """
    Rscript -e "rmarkdown::render('${template_file}', \
      output_file = 'k-${n_metaprograms}_metaprograms-summary-report.html',\
      params = list(metaprograms_object_file = '${metaprograms_file}', \
                    metaprograms_metrics_file = '${metrics_file}', \
                    combined_scores_file = '${combined_score_file}', \
                    category_file = '${gene_set_category_file}', \
                    sample_metadata_file = '${params.sample_metadata_file}'))"
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
  publishDir "${params.outdir}/metaprograms/scored-sces/${sample_id}", mode: 'copy'
  input:
    tuple val(unique_id),
          val(sample_id),
          path(sce_file),
          path(annotations_file)
    path combined_scores_file
  output:
    path("${unique_id}_annotated.rds")
  script:
    """
    05-score-and-annotate-sce.R \
      --sce_file ${sce_file} \
      --unique_id ${unique_id} \
      --annotations_file ${annotations_file} \
      --scores_file ${combined_scores_file} \
      --output_file "${unique_id}_annotated.rds"
    """
  stub:
    """
    touch "${unique_id}_annotated.rds"
    """
}



workflow metaprograms {
  take:
    cnmf_ch  // [sample_id, unique_id, cnmf output, subset_id]
    sce_ch   // [unique_id, sample_id, sce_file, annotations_file]
  main:

    // prepare strings for nanming the metaprogram output directory
    cnmf_k_range_str = params.cnmf_k_range.replaceAll(",", "-")
    filter_str = params.filter_spectra ? "filtered" : "unfiltered"

    // channel of n_metaprograms
    n_metaprograms_ch = channel.fromList(params.n_metaprograms.toString().tokenize(","))

    // channel with all combinations of subset Id and n_metaprograms 
    metaprograms_ch = cnmf_ch
      .groupTuple(by: 3) // group by subset id
       // add in n_metaprograms for each subset id
      .combine(n_metaprograms_ch)
      .map{ sample_id_tuple, unique_id_tuple, cnmf_output_dirs, subset_id, n_metaprograms -> 
        
        // define the publish directory and expected output files
        def metaprograms_publish_dir = "${params.outdir}/metaprograms/${subset_id}/metaprograms-${n_metaprograms}_cnmf-${cnmf_k_range_str}_${filter_str}"
        def mp_file = file("${metaprograms_publish_dir}/k-${n_metaprograms}_metaprograms.rds")
        def combined_scores_file = file("${metaprograms_publish_dir}/k-${n_metaprograms}_metaprogram_scores.tsv.gz")
        def combined_background_scores_file = file("${metaprograms_publish_dir}/k-${n_metaprograms}_background_scores.tsv.gz")
        def metrics_file = file("${metaprograms_publish_dir}/k-${n_metaprograms}_metaprograms-metrics.rds")

        [
          subset_id, 
          n_metaprograms,
          unique_id_tuple,
          cnmf_output_dirs, 
          metaprograms_publish_dir, 
          mp_file, 
          combined_scores_file, 
          combined_background_scores_file, 
          metrics_file
        ]
      }
      // if all the outputs exist for a given n_metaprogram and subset id, let's skip it 
      // because the publishing steps only happen for the full cohort mode, this will only skip for the full cohort 
      .branch { subset_id, n_metaprograms, unique_id_tuple, cnmf_output_dirs, metaprograms_publish_dir, mp_file, combined_scores_file, combined_background_scores_file, metrics_file ->
        cached_metaprograms: 
          !params.repeat_creating_metaprograms &&
          mp_file.exists() && 
          combined_scores_file.exists() && 
          combined_background_scores_file.exists() && 
          metrics_file.exists()
        create_metaprograms: true
      }

    ///////////////////////////////////////////////////////////////////////////////
    //           Metaprogram creation, scoring, and metrics calculation          //
    ///////////////////////////////////////////////////////////////////////////////

    // input for creating metaprograms 
    create_metaprograms_ch = metaprograms_ch.create_metaprograms
      .map{ subset_id, n_metaprograms, unique_id_tuple, cnmf_output_dirs, metaprograms_publish_dir, mp_file, combined_scores_file, combined_background_scores_file, metrics_file -> 
        [subset_id, n_metaprograms, unique_id_tuple, cnmf_output_dirs, metaprograms_publish_dir]
      }

    // generate metaprograms
    generate_metaprograms(create_metaprograms_ch)

    // score metaprograms for each sample
    // make sure to score for each sample for each value of n_metaprograms
    score_input_ch = generate_metaprograms.out
      .transpose(by: 2) // transpose unique id and subset/n metaprogram information 
       // pull out unique id as first element for joining
      .map{ subset_id, n_metaprograms, unique_id, metaprograms_publish_dir, mp_file -> 
        [unique_id, subset_id, n_metaprograms, metaprograms_publish_dir, mp_file]
      }
      // combine with sce channel to get sce file and annotation file information for each unique id in each subset
      // this is why we had to set unique id first in the sce channel
      .combine(sce_ch, by: 0) // combine on unique id

    score_metaprograms(score_input_ch)

    grouped_score_ch = score_metaprograms.out
      .map{ subset_id, n_metaprograms, metaprograms_publish_dir, unique_id, sample_id, score_file, background_score_file -> 
        [subset_id, n_metaprograms, metaprograms_publish_dir, score_file, background_score_file]
      }
      // group by n_metaprograms to get list of score files and background score files for each n_metaprograms
      .groupTuple(by: [0,1,2]) // group by subset_id, n_metaprograms and publish dir


    // create a single data frame for all observed and all background scores
    combine_score_files(grouped_score_ch)

    // calculate metrics 
    metaprogram_metrics(
      generate_metaprograms.out, 
      file(params.custom_genesets_file),
      file(params.msigdb_genesets_file),
      file(params.gene_set_category_file)
    )

    /////////////////////////////////////////////////////
    //               Mix with cached outputs          //
    /////////////////////////////////////////////////////

    
    // construct channels for skipped processes and mix with output from actual processes
    all_metaprograms_ch = metaprograms_ch.cached_metaprograms
      .map{ subset_id, n_metaprograms, unique_id_tuple, cnmf_output_dirs, metaprograms_publish_dir, mp_file, combined_scores_file, combined_background_scores_file, metrics_file ->
        [subset_id, n_metaprograms, unique_id_tuple, metaprograms_publish_dir, mp_file]
      }
      .mix(generate_metaprograms.out)

    all_scores_ch = metaprograms_ch.cached_metaprograms
      .map{ subset_id, n_metaprograms, unique_id_tuple, cnmf_output_dirs, metaprograms_publish_dir, mp_file, combined_scores_file, combined_background_scores_file, metrics_file ->
        [subset_id, n_metaprograms, combined_scores_file, combined_background_scores_file]
      }
      .mix(combine_score_files.out)

    // collect all metrics files for the combined report
    all_metrics_ch = metaprograms_ch.cached_metaprograms
      // match the output of the metaprograms metrics process so we can combine them 
      .map{ subset_id, n_metaprograms, unique_id_tuple, cnmf_output_dirs, metaprograms_publish_dir, mp_file, combined_scores_file, combined_background_scores_file, metrics_file ->
        [subset_id, n_metaprograms, metaprograms_publish_dir, metrics_file]
      }
      .mix(metaprogram_metrics.out)


    /////////////////////////////////////////////////////
    //               Combined metrics report           //
    /////////////////////////////////////////////////////

    // define path to template reports
    individual_report_template_file = "${projectDir}/templates/individual-metaprogram-evaluation.Rmd" 
    combined_report_template_file = "${projectDir}/templates/combined-metaprogram-metrics.Rmd"

    // use the grouped channel as input for the combined and indivudal reports
    grouped_metrics_ch = all_metrics_ch
      // group by subset id to get list of metrics files for each
      .groupTuple(by: 0)
      // drop publish directory since we don't need it 
      .map{ subset_id, n_metaprograms, metaprograms_publish_dir, metrics_files -> 
        [subset_id, n_metaprograms, metrics_files]
      }

    // create input channel for combined report 
    metrics_report_ch = grouped_metrics_ch
      // group scores channel by subset id 
      // then combine with metrics files using subset id
      .combine(all_scores_ch.groupTuple(by:0), by: [0, 1])

    // render combined report
    combined_metaprogram_report(
      metrics_report_ch,
      file(combined_report_template_file)
    )

    /////////////////////////////////////////////////////
    //               Identify optimal k               //
    /////////////////////////////////////////////////////

    // make a channel of sce file, annotation file and score file 
    // only for the optimal k value 
    optimal_k_scores_ch = all_scores_ch
      // drop publish directory and background score files 
      .map{ subset_id, n_metaprograms, combined_score_file, combined_background_score_file -> 
        [subset_id, n_metaprograms, combined_score_file]
      }
      .combine(combined_metaprogram_report.out, by: 0) // combine by subset id to get optimal k for each subset
      // filter to only get score outputs for optimal k value
      .filter{ subset_id, n_metaprograms, combined_score_file, optimal_k, combined_report -> 
        n_metaprograms == optimal_k 
      }
      // drop extra k and report
      .map{ subset_id, n_metaprograms, combined_score_file, optimal_k, combined_report -> 
        [subset_id, optimal_k, combined_score_file]
      }

    // only score and annotate if specified and subset mode is not used
    if (params.score_sces && !params.subset_mode) {
      // pull out just the combined scores file to use for annotating sces 
      combined_scores_file = optimal_k_scores_ch
        .collect{ subset_id, optimal_k, combined_score_file -> combined_score_file }

      // output sces with annotations and scores for the optimal k value
      score_and_annotate_sce(sce_ch, combined_scores_file)
    }

    // only render a report for the optimal k value
    // need a tuple of [subset_id, n_metaprograms, metaprograms_file, metrics_file, combined_score_files] for the report
    report_ch = all_metaprograms_ch
      .combine(combined_metaprogram_report.out, by: 0 ) // combine by subset id and get optimal k for each subset
      .filter{ subset_id, n_metaprograms, unique_id_tuple, metaprograms_publish_dir, metaprograms_file, optimal_k, combined_report -> 
        n_metaprograms == optimal_k 
      }
      .map { subset_id, n_metaprograms, unique_id_tuple, metaprograms_publish_dir, metaprograms_file, optimal_k, combined_report -> 
        [subset_id, optimal_k, metaprograms_file]
      }
      .join(grouped_metrics_ch.transpose(by: [1, 2]), by: [0, 1]) // join on subset_id and n_metaprograms to keep only the optimal k metrics file for each subset
      .join(optimal_k_scores_ch, by: [0,1]) // join on subset_id and optimal_k


    individual_metaprogram_report(
      report_ch,
      file(individual_report_template_file),
      file(params.gene_set_category_file)
    )

  emit:
    metaprograms = all_metaprograms_ch // [subset_id, n_metaprograms, unique_id_tuple, metaprograms_publish_dir, metaprograms_file]
    scores = all_scores_ch // [subset_id, n_metaprograms, combined_score_file, combined_background_score_file]
    metrics = all_metrics_ch // [subset_id, n_metaprograms, metaprograms_publish_dir, metrics_file]
}
