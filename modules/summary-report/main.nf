#!/usr/bin/env nextflow

// Module for generating the summary report 

process render_report {
  container params.ews_nf_container
  tag "${unique_id}"
  label 'mem_8'
  time '2h'
  publishDir "${params.outdir}/reports/${sample_id}"
  input:
    tuple val(unique_id),
          val(sample_id),
          path(sce_file), 
          path(annotations_file)
    path template_file
    path plotting_functions_file
    path validation_groups_file
    path marker_genes_file
    path color_palette_file
  output:
    tuple val(unique_id),
          val(sample_id),
          path(summary_report_html)
  script:
    summary_report_html = "${unique_id}_summary-report.html"
    """
    Rscript -e "rmarkdown::render('${template_file}', \
      output_file = '${summary_report_html}',\
      params = list(unique_id = '${unique_id}', \
                    sce_file = '${sce_file}', \
                    annotations_results_file = '${annotations_file}', \
                    validation_groups_file = '${validation_groups_file}', \
                    marker_genes_file = '${marker_genes_file}', \
                    color_palette_file = '${color_palette_file}', \
                    plotting_functions_file = '${plotting_functions_file}',\
                    marker_gene_cutoff = ${params.marker_gene_cutoff}))"
    """
  stub:
    summary_report_html = "${unique_id}_summary-report.html"
    """
    touch ${summary_report_html}
    """
}


workflow summary_report {
  take:
    library_ch  // [unique_id, sample_id, path to processed SCE file]
    annotations_ch // [unique_id, sample_id, annotations tsv]
    infercnv_ch // [unique id, sample_id, path to cnv heatmap, path to cnv metadata]
  main:

    // define path to template report 
    report_template_file = "${projectDir}/templates/summary-report.Rmd" 
    plotting_functions_file = "${projectDir}/templates/utils/plotting-functions.Rmd"

    // join sce file with aucell results 
    report_ch = library_ch
      .join(annotations_ch, by: [0,1])
      .join(infercnv_ch, by:  [0,1])
      // drop the cnv metadata after joining
      .map{ unique_id, sample_id, sce, annotations_file, cnv_metadata -> tuple(
        unique_id,
        sample_id,
        sce, 
        annotations_file
      )}

    // render report 
    render_report(
      report_ch, 
      file(report_template_file),
      file(plotting_functions_file),
      file(params.validation_groups_file),
      file(params.marker_genes_file),
      file(params.color_palette_file)
    )

  emit:
    render_report.out // [unique_id, sample_id, report file, annotations tsv]
}
