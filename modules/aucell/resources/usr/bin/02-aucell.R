#!/usr/bin/env Rscript

# This script is used to run `AUCell` on a single SCE object for a set of marker gene sets
# gene sets used are custom gene sets and a set of Ewing specific gene sets from MsigDB
# gene sets should be provided as a GeneSetCollection object (as created by 01-build-genesetcollection.R)

# the results are exported as a single TSV file with the following columns:
# `gene_set`, `barcodes`, `auc`, `auc_threshold`, and `unique_id`

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    default = "",
    help = "Path to RDS file containing a processed SingleCellExperiment object to use with AUCell."
  ),
  make_option(
    opt_str = c("--gene_set_collection_file"),
    type = "character",
    default = "",
    help = "Path to RDS file containing GeneSetCollection object to use with AUCell."
  ),
  make_option(
    opt_str = c("--max_rank_threshold"),
    type = "double",
    default = 0.01,
    help = "Proportion of genes detected to set as the `aucMaxRank`. 
      Must be a number between 0 and 1."
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    help = "Path to TSV file where results will be saved"
  ),
  make_option(
    opt_str = c("--is_reference"),
    type = "store_true",
    default = FALSE,
    help = "Specify that the input SCE object contains a collection of reference cells.
      If used, all cells will be labeled as 'reference' in the 'ref_or_query' column of the output 
      and all cells will be used to get the mean and standard deviation for calculating the z-score"
  ), 
  make_option(
    opt_str = c("--ref_aucell_results"),
    type = "character",
    default = "",
    help = "Path to TSV file containing AUCell results for all reference cells. 
      This is required if --is_reference is FALSE and all cells in this file will be used to calculate z-scores."
  ),
  make_option(
    opt_str = c("-t", "--threads"),
    type = "integer",
    default = 4,
    help = "Number of multiprocessing threads to use."
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = 2025,
    help = "A random seed for reproducibility."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# set the seed
set.seed(opt$seed)

# make sure input files exist
stopifnot(
  "sce file does not exist" = file.exists(opt$sce_file),
  "GeneSetCollection file does not exist" = file.exists(opt$gene_set_collection_file),
  "max_rank_threshold must be between 0 and 1" = (opt$max_rank_threshold <= 1 & opt$max_rank_threshold > 0),
  "Output file must end in .tsv" = stringr::str_detect(opt$output_file, "\\.tsv$")
)

# load SCE
suppressPackageStartupMessages({
  library(SingleCellExperiment)
})


# set up multiprocessing params
if (opt$threads > 1) {
  bp_param <- BiocParallel::MulticoreParam(opt$threads)
} else {
  bp_param <- BiocParallel::SerialParam()
}

# make sure directory exists for writing output
output_dir <- dirname(opt$output_file)
fs::dir_create(output_dir)

# read in SCE
sce <- readr::read_rds(opt$sce_file)

# set up the metadata to add to the output TSV 
# based on if the input is a single library or if it's a collection of reference cells
# if its reference cells, then cells from multiple libraries will be present and they should all be labeled "reference"
if(opt$is_reference){
  # save coldata for joining with results later 
  coldata_df <- colData(sce) |> 
    as.data.frame() |> 
    dplyr::select(barcodes, sample_id, unique_id) |> 
    # add cell id for joining and label everything as reference in this case 
    dplyr::mutate(cell_id = glue::glue("{unique_id}_{barcodes}"),
                  ref_or_query = "reference")
  
} else {
  # if the input is a single library, then we run AUCell on that library 
  # this also includes calculating z-scores, which requires AUCell results for a reference group of cells
  
  # check that reference tsv is provided 
  stopifnot("Missing reference_cells_file" = file.exists(opt$reference_cells_file))
  
  # read in ref cells file with aucell results for all ref cells 
  ref_aucell_df <- readr::read_tsv(opt$reference_cells_file)
  
  # define the list of reference cells
  ref_cells <- ref_aucell_df$cell_id |> 
    unique()
  
  # extract unique id to save to results file from sce file name 
  unique_id <- stringr::word(basename(opt$sce_file), sep = "_", 1)
  
  # get the coldata to join with results later
  coldata_df <- data.frame(
    barcodes = colnames(sce),
    sample_id = metadata(sce)$sample_id,
    unique_id = unique_id,
    cell_id = glue::glue("{unique_id}_{barcodes}")
  ) |> 
    # add ref or query column 
    # label based on if cells are part of the reference or query cells
    dplyr::mutate(ref_or_query = dplyr::if_else(cell_id %in% ref_cells, "reference", "query"))
  
  # set names of sce to be cell_id to keep labeling consistent across run modes
  colnames(sce) <- coldata_df$cell_id
}

# remove genes that are not detected from SCE object
genes_to_keep <- rowData(sce)$detected > 0
filtered_sce <- sce[genes_to_keep, ]

# read in gene set collection to use with AUcell
gene_set_collection <- readr::read_rds(opt$gene_set_collection_file)

# Run AUCell -------------------------------------------------------------------

# extract counts matrix
counts_mtx <- counts(filtered_sce)

# check intersection with gene sets
genes_list <- GSEABase::geneIds(gene_set_collection)
overlap_pct <- genes_list |>
  purrr::map_dbl(\(list){
    num_genes <- length(list)
    intersect(rownames(counts_mtx), list) |>
      length() / num_genes
  })

# if any gene sets don't have enough overlap (cutoff is 20%)
# print a message and quit
if (any(overlap_pct <= 0.20)) {
  message("Gene sets do not have at least 20% of genes present in SCE.
          AUCell will not be run.")
  # make empty data frame and save to output file
  data.frame(
    barcodes = colnames(sce),
    gene_set = NA,
    auc = NA,
    auc_thresholds = NA
  ) |>
    readr::write_tsv(opt$output_file)

  # don't run the rest
  quit(save = "no")
}

# get max rank 
max_rank <- ceiling(opt$max_rank_threshold*nrow(counts_mtx))

# run aucell
auc_results <- AUCell::AUCell_run(
  counts_mtx,
  gene_set_collection,
  aucMaxRank = max_rank,
  BPPARAM = bp_param
)

# Get threshold ----------------------------------------------------------------

# get auc threshold for each geneset
auc_thresholds <- AUCell::AUCell_exploreThresholds(
  auc_results,
  assign = TRUE,
  plotHist = FALSE
) |>
  # extract select auc threshold
  purrr::map_dbl(\(results){
    results$aucThr$selected
  })

# put into a data frame for easy joining with all auc values
threshold_df <- data.frame(
  gene_set = names(auc_thresholds),
  auc_threshold = auc_thresholds
)

# Combine and export results ---------------------------------------------------

# create data frame with auc for each cell and each geneset
auc_df <- auc_results@assays@data$AUC |>
  as.data.frame() |>
  tibble::rownames_to_column("gene_set") |>
  tidyr::pivot_longer(!"gene_set",
                      names_to = "cell_id",
                      values_to = "auc"
  ) |>
  # add in threshold column
  dplyr::left_join(threshold_df, by = "gene_set") |>
  # add ids from coldata 
  dplyr::left_join(coldata_df, by = c("cell_id")) |> 
  dplyr::relocate(gene_set, .after = cell_id)

# get all gene set names 
gene_set_names <- unique(auc_df$gene_set)

# calculate z-scores 
# calculate ref mean and sd and get z-score for all cells in each gene set 
zscore_df <- gene_set_names |> 
  purrr::map(\(name){
    
    # subset results to only include auc for the gene set
    gene_set_df <- auc_df |> 
      dplyr::filter(gene_set == name) 

    # if the object is a single library, grab the reference cells from the ref object for the given gene set
    if(!opt$is_reference){
      ref_auc <- ref_aucell_df |> 
        dplyr::filter(gene_set == name) |> 
        dplyr::pull(auc)
    } else {
      # if the object is the reference, then use the auc_df to get all reference cells
      ref_auc <- gene_set_df$auc
    }
    
    # calculate mean and sd based on reference cells
    ref_mean <- mean(ref_auc)
    ref_sd <- sd(ref_auc)
    
    # calculate z-score for each cell 
    gene_set_df <- gene_set_df |> 
      dplyr::mutate(
        auc_zscore = (auc - ref_mean)/ref_sd
      )
    
    return(gene_set_df)
  }) |> 
  dplyr::bind_rows()


# export results as table
readr::write_tsv(zscore_df, opt$output_file)
