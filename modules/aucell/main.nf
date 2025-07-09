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
  tag "${unique_id}"
  label 'mem_24'
  cpus 2
  time '2h'
  publishDir "${params.outdir}/aucell/${sample_id}", mode: 'copy'
  input:
    tuple val(unique_id),
          val(sample_id),
          path(sce_file)
    path geneset_collection_file
    val aucell_max_rank_threshold
  output:
    tuple val(unique_id),
          val(sample_id),
          path(aucell_output_file)
  script:
    aucell_output_file = "${unique_id}_aucell-results.tsv"
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
    aucell_output_file = "${unique_id}_aucell-results.tsv"
    """
    touch ${aucell_output_file}
    """
}

process aucell_zscore {
  container params.ews_nf_container
  tag "aucell-zscore"
  label 'mem_24'
  cpus 2
  time '2h'
  publishDir "${params.outdir}/aucell", mode: 'copy'
  input:
    path results_files
    path ref_annotations
  output:
    path "**/*_aucell-zscores.tsv", emit: query_zscores
    path ref_zscores, emit: ref_zscores 
  script:
    def aucell_results_files = results_files.join(',')
    ref_zscores = "aucell-reference-zscores.tsv"
    """
    03-calculate-zscore.R \
      --aucell_results_files ${aucell_results_files} \
      --reference_cells_file ${ref_annotations} \
      --reference_output_file ${ref_zscores}
    """

  stub:
    ref_zscores = "aucell-reference-zscores.tsv"
    """
    # use stub sample/library ids
    mkdir SCPCS999999
    touch SCPCS999999/SCPCL999999_aucell-zscores.tsv
    touch ${ref_zscores}
    """
}


workflow aucell {
  take:
    library_ch  // [unique_id, sample_id, path to processed SCE file]
    ref_ch // [ref_sce, ref_annotations]
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

    // extract ref annotations from channel to use for zscores
    ref_annotations = ref_ch.map{ ref_sce, ref_annotations -> ref_annotations}

    // grab aucell output files
    zscore_ch = run_aucell.out.map{ sample_id, unique_id, output_file -> output_file }
      .collect()

    // calculate z-scores 
    aucell_zscore(
      zscore_ch,
      ref_annotations
    )

    // pull out individual zscore files to a channel of [unique id, aucell results]
    zscores_ch = aucell_zscore.out.query_zscores
      .flatten()
      .map{ file -> 
        def unique_id = file.name.replace("_aucell-zscores.tsv", "") // pull out unique id from filename
        def sample_id = file.parent.name // sample id is folder name
        tuple(unique_id, sample_id, file)
      }

  emit:
    query = zscores_ch // [unique id, sample_id, aucell zscore output]
    ref = aucell_zscore.out.ref_zscores // [ref_zscore_tsv]
}
