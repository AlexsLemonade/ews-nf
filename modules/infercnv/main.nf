#!/usr/bin/env nextflow

// Module for running inferCNV

process make_gene_order_file {
  container params.ews_nf_container
  label 'mem_8'
  time '2h'
  publishDir "${params.outdir}/reference", mode: 'copy'
  input:
    path gtf_file
    val gene_order_file
  output:
    path gene_order_file
  script:
    gene_order_file = file(gene_order_file).name
    """
    01-make-gene-order-file.R \
      --gtf_file ${gtf_file} \
      --output_file ${gene_order_file}
    """
  stub:
    gene_order_file = file(gene_order_file).name
    """
    touch ${gene_order_file}
    """
}

process run_infercnv {
  container params.ews_nf_container
  tag "${unique_id}"
  label 'big_mem'
  publishDir "${params.outdir}/infercnv/${sample_id}", mode: 'copy'
  input:
    tuple val(unique_id),
          val(sample_id),
          path(sce_file)
    path reference_cells_file
    path gene_order_file
  output:
    tuple val(unique_id),
          val(sample_id),
          path(infercnv_metadata)
  script:
    infercnv_metadata = "${unique_id}_infercnv-metadata.tsv"
    """
    02-run-infercnv.R \
      --sce_file ${sce_file} \
      --unique_id ${unique_id} \
      --reference_file ${reference_cells_file} \
      --gene_order_file ${gene_order_file} \
      --output_metadata_file ${infercnv_metadata} \
      --temp_dir \$PWD \
      --threads ${task.cpus} \
      --seed ${params.seed}
    """

  stub:
    infercnv_metadata = "${unique_id}_infercnv-metadata.tsv"
    """
    touch ${infercnv_metadata}
    """
}


workflow infercnv {
  take:
    library_ch  // [unique_id, sample_id, path to processed SCE file]
  main:

    infercnv_ch = library_ch
      .map{ unique_id, sample_id, sce_file -> tuple(
        unique_id,
        sample_id, 
        sce_file,
        // path to metadata file
        file("${params.outdir}/infercnv/${sample_id}/${unique_id}_infercnv-metadata.tsv")
      )}
      .branch{
        // branch samples based on if infercnv exists and skip_infercnv is used 
        skip_infercnv: 
        (
          params.skip_infercnv 
          && it[3].exists()
        )
        run_infercnv: true
      }

    // run infercnv using a normal reference 
    run_infercnv_ch = infercnv_ch.run_infercnv
      .map{ unique_id, sample_id, sce_file, _metadata_file -> tuple(
        unique_id,
        sample_id,
        sce_file
      )}

    run_infercnv(
      run_infercnv_ch,
      file(params.ref_sce_file, checkIfExists: true), 
      file(params.gene_order_file, checkIfExists: true)
    )


    output_ch = infercnv_ch.skip_infercnv
      // drop sce file from channel where infercnv wasn't run
      .map{ unique_id, sample_id, _sce_file, metadata_file -> tuple(
        unique_id,
        sample_id,
        metadata_file
      )}
      // mix with samples that were run
      .mix(run_infercnv.out)


  emit:
    output_ch // [unique_id, sample_id, infercnv metadata]
}
