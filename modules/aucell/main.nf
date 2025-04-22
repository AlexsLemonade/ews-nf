#!/usr/bin/env nextflow

// Module for running AUCell

process create_geneset_collection {
  container params.ews_nf_container
  label 'mem_8'
  time '2h'
  input:
    path custom_genesets_file
    path msigdb_genesets_file
  output:
    path geneset_collection_file
  script:
    geneset_collection_file = "ews-gene-set-collection.rds"
    """
    01-build-genesetcollection.R \
      --custom_genesets_file ${custom_genesets_file} \
      --msigdb_genesets ${msigdb_genesets_file} \
      --output_file ${geneset_collection_file}
    """
  stub:
    geneset_collection_file = "ews-gene-set-collection.rds"
    """
    touch ${geneset_collection_file}
    """
}

process run_aucell {
  container params.ews_nf_container
  tag "${sample_id}-${library_id}"
  label 'mem_24'
  cpus 2
  time '2h'
  publishDir "${params.outdir}/aucell/${sample_id}", mode: 'copy'
  input:
    tuple val(sample_id),
          val(library_id),
          path(sce_file)
    path geneset_collection_file
    val aucell_max_rank_threshold
  output:
    tuple val(sample_id),
          val(library_id),
          path(aucell_output_file)
  script:
    aucell_output_file = "${library_id}_aucell-results.tsv"
    """
    02-aucell.R \
      --sce_file ${sce_file} \
      --gene_set_collection_file ${geneset_collection_file} \
      --max_rank_threshold ${aucell_max_rank_threshold} \
      --output_file ${aucell_output_file} \
      --threads ${task.cpus} \
      --seed ${params.seed}
    """

  stub:
    aucell_output_file = "${library_id}_aucell-results.tsv"
    """
    touch ${aucell_output_file}
    """
}


workflow aucell {
  take:
    library_ch  // [sample id, library_id, path to processed SCE file]
  main:

    // create gene set collection to use with AUCell
    create_geneset_collection(
      file(params.custom_genesets_file),
      file(params.msigdb_genesets_file)
    )

    // run aucell 
    run_aucell(
      library_ch, 
      create_geneset_collection.out, 
      params.aucell_max_rank_threshold
    )

  emit:
    run_aucell.out // [sample_id, project_id, aucell output file]
}
