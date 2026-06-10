#!/usr/bin/env nextflow

// Module for generating the summary report 

process export_annotations {
  container params.ews_nf_container
  tag "${unique_id}"
  label 'mem_8'
  time '2h'
  publishDir "${params.annotations_output_dir}/${sample_id}", mode: 'copy'
  input:
    tuple val(unique_id),
          val(sample_id),
          path(sce_file), 
          path(aucell_file),
          path(infercnv_metadata)
    path ref_aucell_file
    path gene_set_category_file
  output:
    tuple val(unique_id),
          val(sample_id),
          path(output_annotations_file)
  script:
    output_annotations_file = "${unique_id}_annotations.tsv"
    """
    01-annotate-tumor-cells.R \
      --unique_id ${unique_id} \
      --sce_file "${sce_file}" \
      --aucell_results_file "${aucell_file}" \
      --reference_aucell_file "${ref_aucell_file}" \
      --gene_set_category_file "${gene_set_category_file}" \
      --infercnv_results_file "${infercnv_metadata}" \
      --marker_gene_cutoff ${params.marker_gene_cutoff} \
      --output_file "${output_annotations_file}"
    """
  stub:
    output_annotations_file = "${unique_id}_annotations.tsv"
    """
    touch ${output_annotations_file}
    """
}


workflow annotate_tumor_cells {
  take:
    library_ch  // [unique_id, sample_id, path to processed SCE file]
    aucell_query // [unique id, sample_id, path to aucell results]
    aucell_ref // path to file with all reference cell results from aucell 
    infercnv_ch // [unique id, sample_id, path to cnv metadata]
  main:

    // join sce file with aucell results 
    annotation_ch = library_ch
      .join(aucell_query, by: [0,1])
      .join(infercnv_ch, by: [0,1])
      // drop the extra sample id from joining 
      .map{ unique_id, sample_id, sce, aucell_file, cnv_metadata -> tuple(
        unique_id,
        sample_id,
        sce, 
        aucell_file,
        cnv_metadata
      )}

    // render report 
    export_annotations(
      annotation_ch, 
      aucell_ref,
      file(params.gene_set_category_file)
    )

  emit:
    export_annotations.out // [unique_id, sample_id, annotations tsv]
}
