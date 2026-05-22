#!/usr/bin/env nextflow

// **** Included workflows and processes from modules ****
include { aucell } from './modules/aucell'
include { build_ref } from './modules/build-reference'
include { make_gene_order_file } from './modules/infercnv'
include { infercnv } from './modules/infercnv'
include { annotate_tumor_cells } from './modules/annotate-tumor-cells'
include { summary_report } from './modules/summary-report'
include { cnmf } from './modules/cnmf'
include { metaprograms } from './modules/metaprograms'
include { compare_metaprograms } from './modules/compare-metaprograms'

include { parseProcessedFile } from './lib/utils.nf'

// **** Parameter checks ****
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
param_error = false

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

// use to build all reference files 
workflow reference {

  // create a channel of [unique_id, sample id, path to processed SCE file, number of processed cells]
  // this includes all samples
  library_ch = Channel.fromPath("${params.scpca_nf_dir}/**/**_processed.rds", type: 'file')
    .map{ file -> 
      def (unique_id, sample_id, rds_file, num_processed_cells) = parseProcessedFile(file)

      [unique_id, sample_id, rds_file, num_processed_cells] // emit tuple of unique id, sample id, file path, and number of processed cells
    }
    .branch{ _unique_id, _sample_id, rds_file, num_processed_cells -> 
      // only include samples with at least 200 processed cells
      include: (rds_file.size() > 0 && num_processed_cells >= 200) || rds_file.name.startsWith("STUB")
      skip: true
    }

  // send unique ids for all libraries that don't have enough cells to log 
  library_ch.skip
    .subscribe{
      log.error("There are fewer than 200 cells found in the processed object for ${it[0]} and will not be included when building the normal cell reference.")
    }

  // define tuple of output ref files 
  ref_files = tuple(
    params.ref_sce_file,
    params.ref_annotations_file
  )
  // build reference of normal cells
  include_library_ch = library_ch.include
    .map{ unique_id, sample_id, rds_file, num_processed_cells -> 
      [unique_id, sample_id, rds_file] // emit tuple of unique id, sample id, and file path for included samples
    }
  build_ref(include_library_ch, ref_files)

  // create gene order file to use with inferCNV
  make_gene_order_file(
    file(params.gtf_file), 
    params.gene_order_file
  )
}

// default entry point to run main workflow 
workflow annotation {

  sample_ids = params.sample_ids?.tokenize(';, ') ?: []
  run_all = sample_ids[0] == "All"

  // create a channel of [unique_id, sample id, path to processed SCE file]
  library_ch = Channel.fromPath("${params.scpca_nf_dir}/**/**_processed.rds", type: 'file')
    .map{ file -> 
      def (unique_id, sample_id, rds_file, num_processed_cells) = parseProcessedFile(file)
      
      [unique_id, sample_id, rds_file, num_processed_cells] // emit tuple of unique id, sample id, file path, and number of processed cells
    }
    // filter to only specified samples
    // use either unique id or sample id for filtering
    .filter{ run_all || (it[0] in sample_ids) || (it[1] in sample_ids) }
    .branch{ _unique_id, _sample_id, rds_file, num_processed_cells -> 
      // only include samples with at least 200 processed cells
      run_workflow: (rds_file.size() > 0 && num_processed_cells >= 200) || rds_file.name.startsWith("STUB")
      skip_workflow: true
    }

  // send unique ids for all libraries that don't have enough cells to log 
  library_ch.skip_workflow
    .subscribe{
      log.error("There are fewer than 200 cells found in the processed object for ${it[0]}.")
    }

  run_ch = library_ch.run_workflow
    .map{ unique_id, sample_id, rds_file, num_processed_cells -> 
      [unique_id, sample_id, rds_file] // emit tuple of unique id, sample id, and file path for included samples
    }

  // run aucell 
  aucell(run_ch)
  
  // run infercnv 
  infercnv(run_ch)

  // annotate tumor cells
  annotate_tumor_cells(run_ch, aucell.out.query, aucell.out.ref, infercnv.out)

  // generate report 
  summary_report(run_ch, annotate_tumor_cells.out, infercnv.out)
  
}

// workflow for generating metaprograms using NMF across all samples 
workflow nmf_metaprograms {

  // create a channel of [sample id, unique_id, path to processed H5 file]
  sample_ch = Channel.fromPath("${params.scpca_nf_dir}/**/**_processed_rna.h5ad", type: 'file')
    .map{ file -> 
      def (unique_id, sample_id, h5_file, num_processed_cells) = parseProcessedFile(file, "_processed_rna.h5ad")
      
      // note that here we will use sample id first for easier joining with subset ids in this workflow
      [sample_id, unique_id, h5_file, num_processed_cells]
    }
    .branch{ _sample_id, _unique_id, h5_file, num_processed_cells -> 
      run_workflow: (h5_file.size() > 0 && num_processed_cells >= 200) || h5_file.name.startsWith("STUB")
      skip_workflow: true
    }

  // send unique ids for all libraries that don't have enough cells to log 
  sample_ch.skip_workflow
    .subscribe{
      log.error("There are fewer than 200 cells found in the processed object for ${it[1]}.")
    }

  // if a subset tsv file is provided AND subset mode is enabled
  // create the sample grouping channel with subset ids and sample ids based on the provided tsv file
  if( file(params.subset_tsv).exists() && params.subset_mode ) {

    // define subsets to run based on parameter
    subset_ids = params.subset_ids?.tokenize(';, ') ?: []
    run_all_subsets = subset_ids[0] == "All"

    // set subset ids to run 
    sample_group_ch = Channel.fromPath(params.subset_tsv, type: 'file')
      .splitCsv(header: true, sep: '\t')
      .map{ row -> 
        def sample_id = row['scpca_sample_id']
        def subset_id = row['subset_id']
        [sample_id, subset_id]
      }
      // keep all subsets or filter to specified subset ids to run 
      .filter{ it -> run_all_subsets || it[1] in subset_ids }

  } else {
    // if using the full mode and not subset mode then grab all samples with a processed h5 file 
    // provide a subset id of "full_cohort" to be used to label the publishDir
    sample_group_ch = sample_ch.run_workflow
      .map{ sample_id, unique_id, h5_file, num_processed_cells -> 
        [sample_id, "full_cohort"] 
      }

  }

  // collect sample ids to filter which samples will be run through cNMF
  sample_ids = sample_group_ch.collect{ sample_id, subset_id -> sample_id }

  nmf_run_ch = sample_ch.run_workflow
    .map { sample_id, unique_id, h5_file, num_processed_cells ->
      def annotations_file = file("${params.outdir}/annotations/${sample_id}/${unique_id}_annotations.tsv", checkIfExists: true)
      [sample_id, unique_id, h5_file, annotations_file]
    }
    .filter{ it -> it[0] in sample_ids.getVal() } // filter to only samples included in sample group channel

  // run cnmf on each sample and output 
  cnmf(nmf_run_ch)

  metaprograms_input_ch = cnmf.out
    .combine(sample_group_ch, by: 0) // combine on sample id to add subset id information to cnmf output
    .unique() // make sure we only have one entry per sample/unique id combination

  // channel of sce files and annotations to use as input for scoring metaprograms
  sce_ch = nmf_run_ch.map{ sample_id, unique_id, h5_file, annotations_file ->
    def sce_file = file(h5_file.toUriString().replaceAll(/_processed_rna\.h5ad$/, "_processed.rds"), checkIfExists: true)
    // switch to setting unique ID first since when we join this with the metaprograms we want to use unique id and not just sample id
    [unique_id, sample_id, sce_file, annotations_file]
  }

  // generate and score metaprograms across all samples 
  metaprograms(metaprograms_input_ch, sce_ch)

}

workflow shared_metaprograms {

  // get subset Ids from the subset tsv file and pair with metaprogram output files for each subset
  subset_ch = Channel.fromPath(params.subset_tsv, type: 'file')
    .splitCsv(header: true, sep: '\t')
    .flatMap{ row -> 
      def subset_id = row['subset_id']
      // don't check if exists for subset file because we don't want to throw an error
      def subset_metaprograms_file = file("${params.outdir}/metaprograms/${subset_id}/*_metaprograms.rds", glob: true)
      def full_metaprograms_file = file("${params.outdir}/metaprograms/full_cohort/*_metaprograms.rds", glob: true, checkIfExists: true)
      
      // now if the subset file is empty log it and move on
      if (subset_metaprograms_file.isEmpty()) {
        log.error("No metaprogram output files found for subset ${subset_id} in ${params.outdir}/metaprograms/${subset_id}/.")
        return [] // nothing gets emitted to the channel 
      } 

      // make sure we get a tuple for each subset id 
      // only expect one file for subset and one file for full cohort so make sure we don't return a nested tuple 
      return [[subset_id, subset_metaprograms_file.first(), full_metaprograms_file.first()]]
    }
    .unique() // only one file per subset

  compare_metaprograms(subset_ch)

}


// entry workflow
workflow {

  // check that the workflow is valid
  // if not provided, then default to running the annotation workflow
  def valid_workflows = ["reference", "nmf_metaprograms", "annotation", "shared_metaprograms"]
  if (params.workflow != null && !valid_workflows.contains(params.workflow)) {
    error "Invalid workflow '${params.workflow}'. Must be one of: reference, nmf_metaprograms, annotation, shared_metaprograms."
  }

  if (params.workflow == "reference") {
    reference()
  } else if (params.workflow == "nmf_metaprograms") {
    nmf_metaprograms()
  } else if (params.workflow == "shared_metaprograms") {
    shared_metaprograms()
  } else {
    annotation()
  }

}
