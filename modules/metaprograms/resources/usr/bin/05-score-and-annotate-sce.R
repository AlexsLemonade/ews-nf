#!/usr/bin/env Rscript

# This script is used to join the annotations from ews-nf 
# and the scores for a set metaprograms to the coldata of the input sce file 

# Read in arguments ------------------------------------------------------------

library(optparse)
library(SingleCellExperiment)

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    default = "",
    help = "Path to processed SCE to add scores and annotations to"
  ),
  make_option(
    opt_str = c("--unique_id"),
    type = "character",
    default = "",
    help = "Unique Id for library"
  ),
  make_option(
    opt_str = c("--annotations_file"),
    type = "character",
    default = "",
    help = "Path to annotations file output from ews-nf with tumor cell annotations"
  ),
  make_option(
    opt_str = c("--scores_file"),
    type = "character",
    default = "",
    help = "Path to TSV file containing scores for all metaprograms and all cells across all libraries"
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    default = "annotated.rds", 
    help = "Path to RDS file to save modified SCE"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# make sure inputs are correct
stopifnot(
  "Processed SCE file does not exist" = file.exists(opt$sce_file),
  "Annotations file does not exist" = file.exists(opt$annotations_file),
  "Score file does not exist" = file.exists(opt$scores_file),
  "Output file must end in .rds" = stringr::str_detect(opt$output_file, "\\.rds$")
)

# read in files
sce <- readr::read_rds(opt$sce_file)

# read in annotations and filter to only cells from that library 
annotations_df <- readr::read_tsv(opt$annotations_file) |> 
  dplyr::filter(unique_id == params$unique_id) |> 
  # match the column name for final annotations that is currenlty used in shared data
  dplyr::mutate(ews_nf_celltype_annotation = annotation_tier0) 

# read in scores file
scores_df <- readr::read_tsv(opt$scores_file) |> 
  dplyr::filter(unique_id == params$unique_id)

# expand to have one row per barcode and one column per mp
wide_scores_df <- scores_df |> 
  dplyr::select(-unique_id) |> 
  tidyr::pivot_wider(
    names_from = metaprogram,
    values_from = mp_score,
    names_prefix = "score_"
  )

# join annotations and scores with coldata then add back to object
colData(sce) <- colData(sce) |> 
  as.data.frame() |>
  dplyr::left_join(annotations_df, by = "barcodes") |> 
  dplyr::left_join(wide_scores_df, by = "barcodes") |> 
  DataFrame(row.names = colnames(sce))

readr::write_rds(sce, opt$output_file, compress = "bz2")
