#!/usr/bin/env nextflow

// Module for running inferCNV

process make_gene_order_file {
  container params.ews_nf_container
  label 'mem_8'
  time '2h'
  input:
    path gtf_file
  output:
    path gene_order_file
  script:
    gene_order_file = "${gtf_file.name.replaceFirst(/\.gtf(\.gz)?$/, '')}.gene_order.txt"
    """
    01-make-gene-order-file.R \
      --gtf_file ${gtf_file} \
      --output_file ${gene_order_file}
    """
  stub:
    gene_order_file = "${gtf_file.name.replaceFirst(/\.gtf(\.gz)?$/, '')}.gene_order.txt"
    """
    touch ${gene_order_file}
    """
}

process run_infercnv {
  container params.ews_nf_container
  tag "${sample_id}-${library_id}"
  label 'mem_48'
  cpus 4
  time '6h'
  publishDir "${params.outdir}/infercnv/${sample_id}", mode: 'copy'
  input:
    tuple val(sample_id),
          val(library_id),
          path(sce_file)
    path reference_cells_file
    path gene_order_file
  output:
    tuple val(sample_id),
          val(library_id),
          path(infercnv_metadata),
          path(infercnv_heatmap)
  script:
    infercnv_metadata = "${library_id}_infercnv-metadata.tsv"
    infercnv_heatmap = "${library_id}_infercnv-heatmap.png"
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
    infercnv_metadata = "${library_id}_infercnv-metadata.tsv"
    infercnv_heatmap = "${library_id}_infercnv-heatmap.png"
    """
    touch ${infercnv_metadata}
    touch ${infercnv_heatmap}
    """
}


workflow infercnv {
  take:
    library_ch  // [sample id, library_id, path to processed SCE file]
    ref_ch // [ref_sce, ref_annotations]
  main:

    // extract ref sce from channel 
    ref_sce = ref_ch.map{ ref_sce, ref_annotations -> ref_sce }

    // create gene order file to use with inferCNV
    make_gene_order_file(
      file(params.gtf_file)
    )

    // run infercnv using a normal reference 
    run_infercnv(
      library_ch,
      ref_sce, 
      make_gene_order_file.out
    )

  emit:
    run_infercnv.out // [sample_id, project_id, infercnv metadata, infercnv heatmap]
}
