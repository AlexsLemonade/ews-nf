#!/usr/bin/env Rscript

# script to calculate the background variation and CV for scores from a set of metaprograms
# the gene names in the metaprogram matrix are shuffled 
# and weights for the top 200 genes are used to get scores for each metaprogram in each sample

# The variation and CV are computed for all scores across all samples 
# for a given metaprogram/replicate combination 
# output includes a dataframe with metaprogram, replicate, score_var, and score_cv

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

# shuffle the genes in gene x mp matrix, recalculate the top 200 gene weights
# generate a score for every cell using the quantile normalized counts
permute_scores <- function(mp_mtx, norm_expr, nreps){
  background <- replicate(nreps, {
    # shuffle genes for metaprograms
    rownames(mp_mtx) <- sample(rownames(mp_mtx))
    
    # extract top 200 genes for each mp
    # first turn matrix into a list of vector weights where the names are gene ids
    mp_list <- asplit(mp_mtx, MARGIN = 2)
    # now return the named list of weights, names are genes
    mp_top_list <- mp_list |>
      purrr::map(\(scores){
        scores |>
          sort(decreasing = TRUE) |>
          head(n=200)
      })
    
    # calculate scores
    scores_df <- mp_top_list |>
      purrr::map(\(mps) {score_mp(mps, norm_expr)}) |>
      dplyr::bind_rows(.id = "metaprogram")
    
  }, simplify = FALSE) |> # return as a list not mtx
    dplyr::bind_rows(.id = "replicate")
  
  return(background)
}


# Read in arguments ------------------------------------------------------------

library(optparse)
library(SingleCellExperiment)

option_list <- list(
  make_option(
    opt_str = c("--sce_files"),
    type = "character",
    default = "",
    help = "Comma separated list of processed SCE files"
  ),
  make_option(
    opt_str = c("--annotation_files"),
    type = "character",
    default = "",
    help = "Comma separated list of annotation files"
  ),
  make_option(
    opt_str = c("--metaprograms_file"),
    type = "character",
    default = here::here("ews-nf-out/metaprograms/metaprograms-10_cnmf-5-10-15_filtered/metaprograms.rds"),
    help = "Path to RDS file containing metaprograms"
  ),
  make_option(
    opt_str = c("--nreps"),
    type = "integer",
    default = 1000,
    help = "Number of replicates to use for permutation testing"
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
  "Metaprograms file does not exist" = file.exists(opt$metaprograms_file),
  "Output file must end in .tsv" = stringr::str_detect(opt$output_file, "\\.tsv$")
)

# get a list of all processed objects and pull out unique ids from each file name
sce_files <- stringr::str_split_1(opt$sce_files, ",")
sce_unique_ids <- stringr::word(basename(sce_files), 1, sep = "_")
names(sce_files) <- sce_unique_ids

# get a list of all annotation files and pull out unique ids from each file name
annotation_files <- stringr::str_split_1(opt$annotation_files, ",")
annotation_unique_ids <- stringr::word(basename(annotation_files), 1, sep = "_")
names(annotation_files) <- annotation_unique_ids

# keep only unique IDs that have both processed object and annotations
unique_ids <- intersect(sce_unique_ids, annotation_unique_ids)
# make sure files are in the same order
sce_files <- sce_files[unique_ids]
annotation_files <- annotation_files[unique_ids]

# check that files exist 
stopifnot(
  "Some or all SCE files do not exist" = all(file.exists(sce_files)),
  "Some or all annotations files do not exist" = all(file.exists(annotation_files))
)

# read in the metaprogram object and extract gene universe and mp mtx
mp_object <- readr::read_rds(opt$metaprograms_file)
gene_universe <- mp_object$gene_universe
mp_mtx <- do.call(cbind, mp_object$metaprogram_list)

# Calculate scores and export --------------------------------------------------

# for each sample, shuffle the mp mtx and re-calculate scores
# do this for each sample for nreps 
# then get the variance and CV for each metaprogram and each replicate 
background_df <- purrr::map2(sce_files, annotation_files, 
                             \(sce_file, annotation_file){
                               
                               # read in files
                               sce <- readr::read_rds(sce_file)
                               annotations_df <- readr::read_tsv(annotation_file)
                               
                               # quantile normalize expr matrix 
                               norm_expr <- quantile_norm_tumor_cells(sce, annotations_df, gene_universe)
                               
                               # if there are tumor cells, permute the scores 
                               # shuffle the gene names in the mp matrix, grab top 200 genes and weights
                               # calculate scores 
                               if(ncol(norm_expr) > 0){
                                 background <- permute_scores(mp_mtx, norm_expr, opt$nreps)
                               } else {
                                 background <- data.frame(
                                   barcodes = NA,
                                   mp_score = NA
                                 )
                               }
                               
                               
                               
                             }) |> 
  dplyr::bind_rows(.id = "unique_id") |> 
  # summary stats for each mp and each rep
  dplyr::group_by(metaprogram, replicate) |> 
  dplyr::summarize(
    score_var = var(mp_score, na.rm = TRUE),
    score_cv = sqrt(score_var) / mean(mp_score)
  )

# export output data frame
readr::write_tsv(background_df, opt$output_file)

