#!/usr/bin/env nextflow

// Module for building normal references to use with other modules 

process save_normal_ref {
  container params.ews_nf_container
  // TODO: Remove this once adding in the next process 
  publishDir "${params.outdir}/build-normal-reference/${sample_id}"
  label 'mem_8'
  time '2h'
  input:
  tuple val(sample_id),
        val(library_id),
        path(sce_file)
    path ref_celltypes_file
  output:
    path(normal_ref_rds)
  script:
    normal_ref_rds = "${library_id}_normal-ref.rds"
    """
    01-save-reference-cells.R \
      --sce_file ${sce_file} \
      --reference_celltypes_file ${ref_celltypes_file} \
      --output_file ${normal_ref_rds}
    """
  stub:
    normal_ref_rds = "${library_id}_normal-ref.rds"
    """
    touch ${normal_ref_rds}
    """
}


workflow build_normal_ref {
  take:
    library_ch  // [sample id, library_id, path to processed SCE file]
  main:

    // save ref rds 
    save_normal_ref(
      library_ch, 
      file(params.ref_celltypes_file)
    )

  // TODO: populate with normal reference output from combining step
  //emit:

}
