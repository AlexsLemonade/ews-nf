#!/usr/bin/env nextflow

// Module for generating the summary report 

process render_report {
  container params.ews_nf_container
  label 'mem_8'
  time '2h'
  publishDir "${params.outdir}/reports/${sample_id}"
  input:
    tuple val(unique_id),
          val(sample_id),
          path(sce_file), 
          path(aucell_file),
          path(infercnv_metadata),
          path(infercnv_heatmap)
    path ref_aucell_file
    path ref_annotations_file
    path template_file
    path plotting_functions_file
    path validation_groups_file
    path marker_genes_file
    path color_palette_file
    path gene_set_category_file
  output:
    tuple val(unique_id),
          val(sample_id),
          path(summary_report_html),
          path(output_annotations_file)
  script:
    summary_report_html = "${unique_id}_summary-report.html"
    output_annotations_file = "${unique_id}_annotations.tsv"
    """
    Rscript -e "rmarkdown::render('${template_file}', \
      output_file = '${summary_report_html}',\
      params = list(unique_id = '${unique_id}', \
                    sce_file = '${sce_file}', \
                    validation_groups_file = '${validation_groups_file}', \
                    marker_genes_file = '${marker_genes_file}', \
                    color_palette_file = '${color_palette_file}', \
                    aucell_results_file = '${aucell_file}', \
                    ref_aucell_file = '${ref_aucell_file}', \
                    ref_annotations_file = '${ref_annotations_file}', \
                    gene_set_category_file = '${gene_set_category_file}', \
                    infercnv_heatmap_file = '${infercnv_heatmap}', \
                    infercnv_results_file = '${infercnv_metadata}', \
                    plotting_functions_file = '${plotting_functions_file}',\
                    output_annotations_file = '${output_annotations_file}'))"
    """
  stub:
    summary_report_html = "${unique_id}_summary-report.html"
    output_annotations_file = "${unique_id}_annotations.tsv"
    """
    touch ${summary_report_html}
    touch ${output_annotations_file}
    """
}


workflow summary_report {
  take:
    library_ch  // [unique_id, sample_id, path to processed SCE file]
    aucell_query // [unique id, sample_id, path to aucell results]
    aucell_ref // path to file with all reference cell results from aucell 
    infercnv_ch // [unique id, sample_id, path to cnv heatmap, path to cnv metadata]
  main:

    // define path to template report 
    report_template_file = "${projectDir}/templates/summary-report.Rmd" 
    plotting_functions_file = "${projectDir}/templates/utils/plotting-functions.Rmd"

    // join sce file with aucell results 
    report_ch = library_ch
      .join(aucell_query, by: 0)
      .join(infercnv_ch, by: 0)
      // drop the extra sample id from joining 
      .map{ unique_id, sample_id, sce, sample_id2, aucell_file, sample_id3, cnv_metadata, cnv_heatmap -> tuple(
        unique_id,
        sample_id,
        sce, 
        aucell_file,
        cnv_metadata,
        cnv_heatmap
      )}

    // render report 
    render_report(
      report_ch, 
      aucell_ref,
      file(params.ref_annotations_file, checkIfExists: true),
      file(report_template_file),
      file(plotting_functions_file),
      file(params.validation_groups_file),
      file(params.marker_genes_file),
      file(params.color_palette_file),
      file(params.gene_set_category_file)
    )

  emit:
    render_report.out // [unique_id, sample_id, report file, annotations tsv]
}
