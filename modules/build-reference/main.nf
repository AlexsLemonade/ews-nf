#!/usr/bin/env nextflow

// Module for building normal references to use with other modules 

process save_normal_ref {
  container params.ews_nf_container
  label 'mem_8'
  time '2h'
  input:
    tuple val(sample_id),
          val(library_id),
          path(sce_file)
    path ref_celltypes_file
  output:
    path ref_rds
  script:
    ref_rds = "${library_id}_ref.rds"
    """
    01-save-reference-cells.R \
      --sce_file ${sce_file} \
      --reference_celltypes_file ${ref_celltypes_file} \
      --output_file ${ref_rds}
    """
  stub:
    ref_rds = "${library_id}_ref.rds"
    """
    touch ${ref_rds}
    """
}

process combine_normal_ref {
  container params.ews_nf_container
  publishDir "${params.outdir}/build-reference"
  // TODO: Bump memory here once we have more samples to work with
  label 'mem_8'
  time '2h'
  input:
    path ref_rds_list
  output:
    tuple path(combined_ref_rds),
          path(ref_annotations)
  script:
    combined_ref_rds = "all-ref-sce.rds"
    ref_annotations = "all-ref-annotations.tsv"
    def ref_files = ref_rds_list.join(',')
    """
    02-combine-reference-cells.R \
      --ref_file_list ${ref_files} \
      --output_sce_file ${combined_ref_rds} \
      --output_annotation_file ${ref_annotations}
    """
  stub:
    combined_ref_rds = "all-ref-sce.rds"
    ref_annotations = "all-ref-annotations.tsv"
    """
    touch ${combined_ref_rds}
    touch ${ref_annotations}
    """
}


workflow build_ref {
  take:
    library_ch  // [sample id, library_id, path to processed SCE file]
  main:

    // save ref rds 
    save_normal_ref(
      library_ch, 
      file(params.ref_celltypes_file)
    )

    // combine into one rds 
    combine_normal_ref(save_normal_ref.out.collect())

  emit:
    combine_normal_ref.out // [ref_sce, ref_annotations]

}
