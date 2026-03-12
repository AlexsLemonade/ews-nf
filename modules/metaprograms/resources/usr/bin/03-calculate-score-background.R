#!/usr/bin/env Rscript

# script to calculate the scores for each shuffled metaprogram in a single library 
# shuffled_mps are found in the mp_object
# output is a data frame with unique_id, barcodes, metaprogram, replicate, and mp_score

# Functions --------------------------------------------------------------------
# quantile normalize the logcounts for tumor cells using all genes in the gene universe
quantile_norm_tumor_cells <- function(sce, annotations_df, gene_universe){
  
  # get a vector of tumor cells for subsetting
  tumor_cells <- annotations_df |> 
    dplyr::filter(stringr::str_detect(annotation_tier0, "tumor")) |> 
    dplyr::pull(barcodes)
  
  # get expression matrix for tumor cells and mp genes of interest 
  expr <- logcounts(sce)[gene_universe,tumor_cells] 
  
  # quantile normalize the scores 
  # preprocessCore requires a regular matrix
  norm_expr <- preprocessCore::normalize.quantiles(as.matrix(expr))
  # make sure to add back cell and gene info
  colnames(norm_expr) <- colnames(expr)
  rownames(norm_expr) <- rownames(expr) 
  
  return(norm_expr)
}

# function for getting a score for a single MP for all cells in an object
score_mp <- function(mp_weights, norm_expr_mtx){
  
  # make sure gene names are present and in the object
  stopifnot(
    "Provided weights are missing gene names" = !is.null(names(mp_weights)),
    "genes are not found in the SCE object" = all(names(mp_weights) %in% rownames(norm_expr_mtx))
  )
  
  # mp genes 
  mp_genes <- names(mp_weights)
  
  # expression subset: genes × cells
  sub_expr <- norm_expr_mtx[mp_genes, ]
  
  # compute dot product for each cell (matrix multiply)
  scores <- as.numeric(mp_weights %*% sub_expr)
  
  df <- data.frame(
    barcodes = colnames(norm_expr_mtx),
    mp_score = scores
  ) 
  
  return(df)
}

# Read in arguments ------------------------------------------------------------

library(optparse)
library(SingleCellExperiment)

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    default = "scpca-nf-out/SCPCP000026/SCPCS000841/SCPCL001228-SCPCS000841_processed.rds",
    help = "processed SCE file"
  ),
  make_option(
    opt_str = c("--annotations_file"),
    type = "character",
    default = "ews-nf-out/annotations/SCPCS000841/SCPCL001228-SCPCS000841_annotations.tsv",
    help = "TSV with tumor cell type annotations"
  ),
  make_option(
    opt_str = c("--metaprograms_file"),
    type = "character",
    default = here::here("ews-nf-out/metaprograms/metaprograms-10_cnmf-5-10-15_filtered/metaprograms.rds"),
    help = "Path to RDS file containing metaprograms"
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    default = "scores-background-stats.tsv", 
    help = "Path to TSV file to save summarized stats calculated from the scores background"
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = 2026,
    help = "A random seed for reproducibility."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# set the seed
set.seed(opt$seed)

# make sure inputs are correct
stopifnot(
  "SCE file does not exist" = file.exists(opt$sce_file),
  "Annotations TSV file does not exist" = file.exists(opt$annotations_file),
  "Metaprograms file does not exist" = file.exists(opt$metaprograms_file),
  "Output file must end in .tsv" = stringr::str_detect(opt$output_file, "\\.tsv$")
)

# read in files
sce <- readr::read_rds(sce_file)
annotations_df <- readr::read_tsv(annotation_file)

# read in the metaprogram object and extract gene universe and shuffled mps
mp_object <- readr::read_rds(opt$metaprograms_file)
gene_universe <- mp_object$gene_universe
shuffled_mp_list <- mp_object$shuffled_metaprograms

# extract unique id from file name
unique_id <- stringr::word(basename(opt$sce_file), 1, sep = "_")

# Calculate scores and export --------------------------------------------------

# quantile normalize gene expression for all tumor cells and all genes in the universe
norm_expr <- quantile_norm_tumor_cells(sce, annotations_df, gene_universe)

# if there are tumor cells, calculate the scores for each permutation 
if(ncol(norm_expr) > 0){
  # for each of the shuffled list of mp genes calculate scores
  background_df <- shuffled_mps |> 
    purrr::map(\(mp_top_list){
      # calculate scores
      scores_df <- mp_top_list |>
        purrr::map(\(mps) {score_mp(mps, norm_expr)}) |>
        dplyr::bind_rows(.id = "metaprogram")
    }) |> 
    dplyr::bind_rows(.id = "replicate") |> 
    dplyr::mutate(unique_id = unique_id)
} else {
  # make an emtpy data frame to be exported 
  background_df <- data.frame(
    barcodes = NA,
    mp_score = NA,
    replicate = NA,
    unique_id = NA
  )
}

# export output data frame
readr::write_tsv(background_df, opt$output_file)

