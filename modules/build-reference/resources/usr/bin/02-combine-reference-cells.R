#!/usr/bin/env Rscript

# This script is used to read in and combine all SCEs containing normal cells into a single SCE to use as a normal reference

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--ref_file_list"),
    type = "character",
    help = "Comma separated list of paths to RDS files containing SingleCellExperiment objects with reference cells only."
  ),
  make_option(
    opt_str = c("--output_sce_file"),
    type = "character",
    help = "Path to RDS file to save combined SingleCellExperiment object with counts from all reference cells."
  ),
  make_option(
    opt_str = c("--output_annotation_file"),
    type = "character",
    help = "Path to TSV file to save table with annotations for all reference cells."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# make sure all ref files exist and output file is an rds file 
stopifnot(
  "Output file must end in .rds" = stringr::str_detect(opt$output_file, "\\.rds$"),
  "List of input SCEs to combine is missing" = !is.null(opt$ref_file_list)
)

# make list of ref files 
ref_sce_files <- unlist(stringr::str_split(opt$ref_file_list, ","))
# check that all files exist 
stopifnot(
  "Not all input SCE files exist" = all(file.exists(ref_sce_files))
)

# load SCE
suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

# Combine SCE and export -------------------------------------------------------

# read in sce 
ref_sce_list <- ref_sce_files |> 
  purrr::map(readr::read_rds) 
# combine into a single reference sce 
combined_ref_sce <- do.call(cbind, ref_sce_list)

# combine annotations into one tsv file and save separately
annotations_df <- colData(combined_ref_sce) |>
  as.data.frame() |> 
  # save column with combined library and cell barcodes
  tibble::rownames_to_column("cell_id")

# export annotations and sce 
readr::write_tsv(annotations_df, opt$output_annotation_file)
readr::write_rds(combined_ref_sce, opt$output_sce_file)

