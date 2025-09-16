#!/usr/bin/env nextflow

// Module for building normal references to use with other modules 

process save_normal_ref {
  container params.ews_nf_container
  label 'mem_8'
  time '2h'
  input:
    tuple val(unique_id),
          val(sample_id),
          path(sce_file)
    path ref_celltypes_file
  output:
    path ref_rds
  script:
    ref_rds = "${unique_id}_ref.rds"
    """
    01-save-reference-cells.R \
      --sce_file ${sce_file} \
      --reference_celltypes_file ${ref_celltypes_file} \
      --output_file ${ref_rds}
    """
  stub:
    ref_rds = "${unique_id}_ref.rds"
    """
    touch ${ref_rds}
    """
}

process combine_normal_ref {
  container params.ews_nf_container
  publishDir "${params.outdir}/reference"
  // TODO: Bump memory here once we have more samples to work with
  label 'mem_8'
  time '2h'
  input:
    path ref_rds_list
    tuple val(ref_sce_file),
          val(ref_annotations_file)
  output:
    tuple path(ref_sce_file),
          path(ref_annotations_file)
  script:
    ref_sce_file = file(ref_sce_file).name
    ref_annotations_file = file(ref_annotations_file).name
    def ref_files = ref_rds_list.join(',')
    """
    02-combine-reference-cells.R \
      --ref_file_list ${ref_files} \
      --output_sce_file ${ref_sce_file} \
      --output_annotation_file ${ref_annotations_file}
    """
  stub:
    ref_sce_file = file(ref_sce_file).name
    ref_annotations_file = file(ref_annotations_file).name
    """
    touch ${ref_sce_file}
    touch ${ref_annotations_file}
    """
}


workflow build_ref {
  take:
    library_ch  // [unique_id, sample_id, path to processed SCE file]
    ref_files // [ref_sce_file, ref_annotations_file]
  main:

    // save ref rds 
    save_normal_ref(
      library_ch, 
      file(params.ref_celltypes_file)
    )

    // combine into one rds 
    combine_normal_ref(save_normal_ref.out.collect(), ref_files)

  emit:
    combine_normal_ref.out // [ref_sce, ref_annotations]

}
