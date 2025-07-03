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
    help = "Path to RDS file containing a processed SingleCellExperiment object to use with AUCell."
  ),
  make_option(
    opt_str = c("--gene_set_collection_file"),
    type = "character",
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

# make sure input files exist
stopifnot(
  "sce file must be specified using `--sce_file`" = !is.null(opt$sce_file),
  "GeneSetCollection file must be specified using `--gene_set_collection_file`" = !is.null(opt$gene_set_collection_file)
)

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

# grab library and sample ID to save to results file in a unique_id column 
library_id <- metadata(sce)$library_id
sample_id <- metadata(sce)$sample_id

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
    names_to = "barcodes",
    values_to = "auc"
  ) |>
  # add in threshold column
  dplyr::left_join(threshold_df, by = "gene_set") |>
  dplyr::relocate(gene_set, .after = barcodes) |>
  dplyr::mutate(
    unique_id = glue::glue("{library_id}-{sample_id}")
  )

# export results as table
readr::write_tsv(auc_df, opt$output_file)
