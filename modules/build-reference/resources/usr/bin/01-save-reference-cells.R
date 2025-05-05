#!/usr/bin/env Rscript

# This script is used to extract any normal cells to be used as a reference 
# for other modules like AUCell and inferCNV
# This includes any cells with an assigned consensus cell type present in the list of reference cell types
# specified with the `--reference_celltypes_file` argument 

# An SCE object with just the counts from reference cells 
# and associated consensus annotations stored in the colData will be saved to an RDS file 

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    help = "Path to RDS file containing a processed SingleCellExperiment object to extract reference cells from."
  ),
  make_option(
    opt_str = c("--reference_celltypes_file"),
    type = "character",
    help = "Path to TSV with list of consensus annotations to include in the reference cells."
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    help = "Path to RDS file to save SingleCellExperiment object with counts from reference cells."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# make sure input files exist
stopifnot(
  "sce file must be specified using `--sce_file`" = !is.null(opt$sce_file),
  "Reference cell type file must be specified using `--reference_celltypes_file`" = !is.null(opt$reference_celltypes_file),
  "Output file must end in .rds" = stringr::str_detect(opt$output_file, "\\.rds$")
)

stopifnot(
  "sce file does not exist" = file.exists(opt$sce_file),
  "Reference cell type file does not exist" = file.exists(opt$reference_celltypes_file)
)

# load SCE
suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

# make sure directory exists for writing output
output_dir <- dirname(opt$output_file)
fs::dir_create(output_dir)

# read in sce file 
sce <- readr::read_rds(opt$sce_file)

# read in reference cell types 
ref_celltypes <- readr::read_tsv(opt$reference_celltypes_file) |> 
  dplyr::pull(consensus_annotation)

# Extract reference cells ------------------------------------------------------

# get library id 
library_id <- metadata(sce)$library_id

# make df with barcodes and consensus cell types for all ref cells 
ref_celltype_df <- colData(sce) |> 
  as.data.frame() |> 
  dplyr::select(barcodes, consensus_celltype_annotation, consensus_celltype_ontology) |> 
  dplyr::mutate(library_id = library_id) |> 
  dplyr::filter(consensus_celltype_annotation %in% ref_celltypes) 

# subset reference cells 
ref_cells <- ref_celltype_df$barcodes
ref_sce <- sce[, ref_cells]

# Build new sce and export -----------------------------------------------------

# reconstruct an SCE with just counts and consensus cell type 
small_ref_sce <- SingleCellExperiment(
  assays = list(counts = counts(ref_sce)), 
  colData = ref_celltype_df
)

# print out a message if no normal cells are found 
if(ncol(small_ref_sce) == 0){
  message("No normal cells identified in provided SingleCellExperiment object")
}

# add library id to colnames 
colnames(small_ref_sce) <- glue::glue("{library_id}-{colnames(small_ref_sce)}")

# save to rds file
readr::write_rds(small_ref_sce, opt$output_file, compress = "gz")
