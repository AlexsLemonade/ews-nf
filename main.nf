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

// **** Parameter checks ****
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
param_error = false

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

// use to build all reference files 
workflow reference {

  // create a channel of [unique_id, sample id, path to processed SCE file]
  // this includes all samples
  library_ch = Channel.fromPath("${params.scpca_nf_dir}/**/**_processed.rds", type: 'file')
    .map{ file -> tuple(
      file.name.replace("_processed.rds", ""), // extract unique id which will be either library_id or <library_id>-<sample_id>
      file.parent.name, // extract sample id from file path
      file // full path to file
    )}
    .branch{
      include: it[2].size() > 0 || it[2].name.startsWith("STUB")
      skip: true
    }

  // send unique ids for all libraries that don't have enough cells to log 
  library_ch.skip
    .subscribe{
      log.error("There are no cells found in the processed object for ${it[0]} and will not be included when building the normal cell reference.")
    }

  // define tuple of output ref files 
  ref_files = tuple(
    params.ref_sce_file,
    params.ref_annotations_file
  )
  // build reference of normal cells
  build_ref(library_ch.include, ref_files)

  // create gene order file to use with inferCNV
  make_gene_order_file(
    file(params.gtf_file), 
    params.gene_order_file
  )
}

// default entry point to run main workflow 
workflow {
  sample_ids = params.sample_ids?.tokenize(';, ') ?: []
  run_all = sample_ids[0] == "All"

  // create a channel of [unique_id, sample id, path to processed SCE file]
  library_ch = Channel.fromPath("${params.scpca_nf_dir}/**/**_processed.rds", type: 'file')
    .map{ file -> tuple(
      file.name.replace("_processed.rds", ""), // extract unique id which will be either library_id or <library_id>-<sample_id>
      file.parent.name, // extract sample id from file path
      file // full path to file
    )}
    // filter to only specified samples
    // use either unique id or sample id for filtering
    .filter{ run_all || (it[0] in sample_ids) || (it[1] in sample_ids) }
    .branch{
      run_workflow: it[2].size() > 0 || it[2].name.startsWith("STUB")
      skip_workflow: true
    }

  // send unique ids for all libraries that don't have enough cells to log 
  library_ch.skip_workflow
    .subscribe{
      log.error("There are no cells found in the processed object for ${it[0]}.")
    }

  run_ch = library_ch.run_workflow

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
  sample_ids = params.sample_ids?.tokenize(';, ') ?: []
  run_all = sample_ids[0] == "All"

  // create a channel of [unique_id, sample id, path to processed H5 file]
  sample_ch = Channel.fromPath("${params.scpca_nf_dir}/**/**_processed_rna.h5ad", type: 'file')
    .map{ file -> tuple(
      file.name.replace("_processed_rna.h5ad", ""), // extract unique id which will be either library_id or <library_id>-<sample_id>
      file.parent.name, // extract sample id from file path
      file // full path to file
    )}
    // filter to only specified samples
    // use either unique id or sample id for filtering
    .filter{ unique_id, sample_id, _h5_file -> 
      run_all || (unique_id in sample_ids) || (sample_id in sample_ids) 
    }
    .branch{ _unique_id, _sample_id, h5_file -> 
      run_workflow: h5_file.size() > 0 || h5_file.name.startsWith("STUB")
      skip_workflow: true
    }

  // send unique ids for all libraries that don't have enough cells to log 
  sample_ch.skip_workflow
    .subscribe{
      log.error("There are no cells found in the processed object for ${it[0]}.")
    }

  nmf_run_ch = sample_ch.run_workflow
    .map { unique_id, sample_id, h5_file ->
      def annotations_file = file("${params.outdir}/annotations/${sample_id}/${unique_id}_annotations.tsv", checkIfExists: true)
      [unique_id, sample_id, h5_file, annotations_file]
    }

  // run cnmf on each sample and output 
  cnmf(nmf_run_ch)

  // generate metaprograms 
  metaprograms(cnmf.out)

}
