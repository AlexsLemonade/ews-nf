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
workflow {
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
  sample_ids = params.sample_ids?.tokenize(';, ') ?: []
  run_all = sample_ids[0] == "All"

  // create a channel of [unique_id, sample id, path to processed H5 file]
  sample_ch = Channel.fromPath("${params.scpca_nf_dir}/**/**_processed_rna.h5ad", type: 'file')
    .map{ file -> 
      def (unique_id, sample_id, h5_file, num_processed_cells) = parseProcessedFile(file, "_processed_rna.h5ad")
      
      [unique_id, sample_id, h5_file, num_processed_cells] // emit tuple of unique id, sample id, file path, and number of processed cells
    }
    // filter to only specified samples
    // use either unique id or sample id for filtering
    .filter{ unique_id, sample_id, _h5_file, _num_processed_cells -> 
      run_all || (unique_id in sample_ids) || (sample_id in sample_ids) 
    }
    .branch{ _unique_id, _sample_id, h5_file, num_processed_cells -> 
      run_workflow: (h5_file.size() > 0 && num_processed_cells >= 200) || h5_file.name.startsWith("STUB")
      skip_workflow: true
    }

  // send unique ids for all libraries that don't have enough cells to log 
  sample_ch.skip_workflow
    .subscribe{
      log.error("There are fewer than 200 cells found in the processed object for ${it[0]}.")
    }

  nmf_run_ch = sample_ch.run_workflow
    .map { unique_id, sample_id, h5_file, num_processed_cells ->
      def annotations_file = file("${params.outdir}/annotations/${sample_id}/${unique_id}_annotations.tsv", checkIfExists: true)
      [unique_id, sample_id, h5_file, annotations_file]
    }

  // run cnmf on each sample and output 
  cnmf(nmf_run_ch)

  // channel of sce files and annotations to use as input for scoring metaprograms
  sce_ch = nmf_run_ch.map{ unique_id, sample_id, h5_file, annotations_file ->
    def sce_file = file(h5_file.toUriString().replaceAll(/_processed_rna\.h5ad$/, "_processed.rds"), checkIfExists: true)
    [unique_id, sample_id, sce_file, annotations_file]
  }

  // generate and score metaprograms across all samples 
  metaprograms(cnmf.out, sce_ch)

}
