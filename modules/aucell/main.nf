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

process build_aucell_reference {
  container params.ews_nf_container
  label 'mem_24'
  cpus 2
  time '2h'
  publishDir "${params.outdir}/aucell", mode: 'copy'
  input:
    path geneset_collection_file
    path ref_sce_file
  output:
    path("aucell-reference-zscores.tsv")
  script:
    """
    02-aucell.R \
      --sce_file ${ref_sce_file} \
      --gene_set_collection_file ${geneset_collection_file} \
      --max_rank_threshold ${params.aucell_max_rank_threshold} \
      --output_file "aucell-reference-zscores.tsv" \
      --is_reference \
      --threads ${task.cpus} \
      --seed ${params.seed}
    """

  stub:
    """
    touch aucell-reference-zscores.tsv
    """
}

process run_aucell {
  container params.ews_nf_container
  tag "${unique_id}"
  label 'mem_24'
  cpus 2
  time '2h'
  publishDir "${params.outdir}/aucell/${sample_id}", mode: 'copy'
  input:
    tuple val(unique_id),
          val(sample_id),
          path(sce_file)
    path ref_aucell_results
    path geneset_collection_file
  output:
    tuple val(unique_id),
          val(sample_id),
          path("*_aucell-results.tsv")
  script:
    """
    02-aucell.R \
      --sce_file ${sce_file} \
      --gene_set_collection_file ${geneset_collection_file} \
      --max_rank_threshold ${params.aucell_max_rank_threshold} \
      --output_file "${unique_id}_aucell-results.tsv" \
      --ref_aucell_results ${ref_aucell_results} \
      --threads ${task.cpus} \
      --seed ${params.seed}
    """

  stub:
    """
    touch "${unique_id}_aucell-results.tsv"
    """
}

workflow aucell {
  take:
    library_ch  // [unique_id, sample_id, path to processed SCE file]
  main:

    // create gene set collection to use with AUCell
    create_geneset_collection(
      file(params.custom_genesets_file),
      file(params.msigdb_genesets_file)
    )

    // build aucell reference file 
    build_aucell_reference(
      create_geneset_collection.out,
      file(params.ref_sce_file)
    )

    // run aucell on each library, including calculating z-scores
    run_aucell(
      library_ch, 
      build_aucell_reference.out,
      create_geneset_collection.out
    )

  emit:
    query = run_aucell.out // [unique id, sample_id, aucell zscore output]
    ref = build_aucell_reference.out // [ref_zscore_tsv]
}
