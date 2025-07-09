#!/usr/bin/env nextflow

// Module for generating the summary report 

process render_report {
  container params.ews_nf_container
  label 'mem_8'
  time '2h'
  publishDir "${params.outdir}/reports/${sample_id}"
  input:
    tuple val(sample_id),
          val(unique_id),
          path(sce_file)
    path template_file
    path validation_groups_file
    path marker_genes_file
    path color_palette_file
  output:
    tuple val(sample_id),
          val(unique_id),
          path(summary_report_html)
  script:
    summary_report_html = "${unique_id}_summary-report.html"
    """
    Rscript -e "rmarkdown::render('${template_file}', \
      output_file = '${summary_report_html}',\
      params = list(unique_id = '${unique_id}', \
                    sce_file = '${sce_file}', \
                    validation_groups_file = '${validation_groups_file}', \
                    marker_genes_file = '${marker_genes_file}', \
                    color_palette_file = '${color_palette_file}'))"
    """
  stub:
    summary_report_html = "${unique_id}_summary-report.html"
    """
    touch ${summary_report_html}
    """
}


workflow summary_report {
  take:
    library_ch  // [sample id, unique_id, path to processed SCE file]
  main:

    // define path to template report 
    report_template_file = "${projectDir}/templates/summary-report.Rmd" 

    // render report 
    render_report(
      library_ch,
      file(report_template_file),
      file(params.validation_groups_file),
      file(params.marker_genes_file),
      file(params.color_palette_file)
    )

  emit:
    render_report.out // [sample_id, unique_id, report file]
}
