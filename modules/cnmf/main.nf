#!/usr/bin/env nextflow

// Module for running cnmf on each sample

process run_cnmf {
  container params.ews_nf_container
  tag "${unique_id}"
  label 'mem_96'
  cpus 4
  time '12h'
  publishDir "${params.outdir}/cnmf/${sample_id}", mode: 'copy'
  input:
    tuple val(unique_id),
          val(sample_id),
          path(h5_file),
          path(annotations_file)
  output:
    tuple val(unique_id),
          val(sample_id),
          path(output_dir)
  script:
    // use the sample ID as the output directory 
    output_dir = sample_id
    """
    01-run-cnmf.py \
      --h5ad_file ${h5_file} \
      --output ${output_dir} \
      --unique_id ${unique_id} \
      --k_components_lower ${params.k_components_lower} \
      --k_components_upper ${params.k_components_upper} \
      --jobs ${task.cpus} \
      --seed ${params.seed}
    """

  stub:
    output_dir = sample_id
    """
    touch ${output_dir}
    """
}


workflow cnmf {
  take:
    sample_ch  // [unique_id, sample_id, path to processed H5 file, annotations file]
  main:

    cnmf_ch = sample_ch
      .map{ unique_id, sample_id, h5_file, annotations_file -> 
        [
          unique_id,
          sample_id, 
          h5_file,
          annotations_file,
          // path to cnmf output
          file("${params.outdir}/cnmf/${sample_id}/${unique_id}", type: "dir")
        ]
      }
      .branch{
        // branch samples based on if cnmf results exist
        skip_cnmf: 
        (
          params.skip_cnmf 
          && it[5].exists() 
        )
        run_cnmf: true
      }

    // run cnmf
    run_cnmf_ch = cnmf_ch.run_cnmf
      .map{ it -> it.dropRight(1) } // drop the cnmf output directory

    run_cnmf(
      run_cnmf_ch
    )

    output_ch = cnmf_ch.skip_cnmf
      // drop sce file from channel where infercnv wasn't run
      .map{ unique_id, sample_id, _h5_file, _annotations_file, cnmf_output_dir -> tuple(
        unique_id,
        sample_id,
        cnmf_output_dir
      )}
      // mix with samples that were run
      .mix(run_cnmf.out)


  emit:
    output_ch // [unique_id, sample_id, cnmf output]
}
