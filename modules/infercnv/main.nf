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
  label 'mem_96'
  cpus 8
  time '6h'
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
          path(infercnv_metadata),
          path(infercnv_heatmap)
  script:
    infercnv_metadata = "${unique_id}_infercnv-metadata.tsv"
    infercnv_heatmap = "${unique_id}_infercnv-heatmap.png"
    """
    02-run-infercnv.R \
      --sce_file ${sce_file} \
      --reference_file ${reference_cells_file} \
      --gene_order_file ${gene_order_file} \
      --output_metadata_file ${infercnv_metadata} \
      --output_heatmap_file ${infercnv_heatmap} \
      --temp_dir \$PWD \
      --threads ${task.cpus} \
      --seed ${params.seed}
    """

  stub:
    infercnv_metadata = "${unique_id}_infercnv-metadata.tsv"
    infercnv_heatmap = "${unique_id}_infercnv-heatmap.png"
    """
    touch ${infercnv_metadata}
    touch ${infercnv_heatmap}
    """
}


workflow infercnv {
  take:
    library_ch  // [unique_id, sample_id, path to processed SCE file]
  main:

    // run infercnv using a normal reference 
    run_infercnv(
      library_ch,
      file(params.ref_sce_file, checkIfExists: true), 
      file(params.gene_order_file, checkIfExists: true)
    )

  emit:
    run_infercnv.out // [unique_id, sample_id, infercnv metadata, infercnv heatmap]
}
