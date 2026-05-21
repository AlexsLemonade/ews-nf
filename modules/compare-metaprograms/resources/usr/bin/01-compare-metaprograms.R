#!/usr/bin/env Rscript

# This script compares all metaprograms in a subset of samples to all metaprograms in the full cohort 
# This includes overall correlation (correlation of all gene weights), 
# correlation of the top 200 gene weights, and the Jaccard index for the top 200 gene ids 

# Much of the metrics calculations were copied from
# https://github.com/AlexsLemonade/ews-nf/blob/main/exploratory-notebooks/10-shared-metaprogram-metrics.Rmd

# Functions --------------------------------------------------------------------

# function to extract the top genes from each MP 
# returns a list of vectors with the gene weights for each MP
# names are Ensembl gene ids
# grabs top 200 genes by default
extract_top_genes <- function(mp_list, num_genes = 200){
  
  # get the top X number of genes for each MP
  mp_top_list <- mp_list |> 
    purrr::map(\(weights){
      weights |> 
        sort(decreasing = TRUE) |> 
        head(n=num_genes)
    })
  
  return(mp_top_list)
}

# helper function for getting the jaccard index between two vectors
calc_jaccard <- function(vec1, vec2){
  jaccard <- length(intersect(vec1, vec2)) / length(union(vec1, vec2))
  return(jaccard)
}

# Read in arguments ------------------------------------------------------------

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--subset_id"),
    type = "character",
    default = "",
    help = "Subset ID"
  ),
  make_option(
    opt_str = c("--subset_metaprograms_file"),
    type = "character",
    default = "",
    help = "Path to the file with metaprograms from a given subset"
  ),
  make_option(
    opt_str = c("--full_metaprograms_file"),
    type = "character",
    default = "",
    help = "Path to the file with metaprograms from the full cohort"
  ),
  make_option(
    opt_str = c("--num_top_genes"),
    type = "integer",
    default = 200,
    help = "Number of top genes to use for calculating correlation and Jaccard index"
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    default = "metaprogram-comparison-metrics.tsv", 
    help = "Path to the tsv file to save the comparison metrics"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# make sure inputs are correct
stopifnot(
  "Subset metaprograms file does not exist" = file.exists(opt$subset_metaprograms_file),
  "Full metaprograms file does not exist" = file.exists(opt$full_metaprograms_file),
  "Output file must end in .tsv" = stringr::str_detect(opt$output_file, "\\.tsv$")
)

# read input files 
subset_mp_obj <- readr::read_rds(opt$subset_metaprograms_file)
full_mp_obj <- readr::read_rds(opt$full_metaprograms_file)

# extract the metaprograms list 
# named vector of gene weights, where names are ensembl ids 
subset_mps <- subset_mp_obj$metaprogram_list
full_mps <- full_mp_obj$metaprogram_list

# get all possible pairs of metaprograms
pairs_df <- expand.grid(
  subset_mp = names(subset_mps),
  full_mp = names(full_mps),
  stringsAsFactors = FALSE
)

# Overall correlation of all weights -------------------------------------------

# correlation of all gene weights
overall_correlation_df <- pairs_df |>
  dplyr::mutate(
    overall_correlation = purrr::map2_dbl(
      subset_mp, full_mp,
      \(x, y) {
        
        # get the actual metaprograms from the names stored in pairs_df
        subset <- subset_mps[[x]]
        full <- full_mps[[y]]
        
        # make sure metaprograms are in the right order and then calculate correlation
        genes <- union(names(subset), names(full))
        subset <- subset[genes]
        full <- full[genes]
        cor(subset, full, method = "pearson", use = "pairwise.complete.obs")
        
      }
    )
  )

# Correlation of top 200 genes -------------------------------------------------

# get top 200 genes and gene weights 
subset_top_genes <- extract_top_genes(subset_mps, opt$num_top_genes)
full_top_genes <- extract_top_genes(full_mps, opt$num_top_genes)

# correlation of just top genes 
top_correlation_df <- pairs_df |>
  dplyr::mutate(
    top_correlation = purrr::map2_dbl(
      subset_mp, full_mp,
      \(x, y) {
        
        # get the actual metaprograms from the names stored in pairs_df
        subset <- subset_top_genes[[x]]
        full <- full_top_genes[[y]]
        
        # make sure metaprograms are in the right order and then calculate correlation
        genes <- union(names(subset), names(full))
        subset <- subset[genes]
        full <- full[genes]
        cor(subset, full, method = "pearson", use = "pairwise.complete.obs")
        
      }
    )
  )

# Jaccard index for top 200 genes ----------------------------------------------

# for each pair, get the jaccard index for gene ids in the top X number of genes
jaccard_df <- pairs_df |>
  dplyr::mutate(
    jaccard = purrr::map2_dbl(
      subset_mp, full_mp,
      \(x, y) {
        
        subset <- names(subset_top_genes[[x]])
        full <- names(full_top_genes[[y]])
        
        calc_jaccard(subset, full)
        
      }
    )
  )

# Combine and export -----------------------------------------------------------

all_stats_df <- overall_correlation_df |> 
  dplyr::left_join(top_correlation_df, by = c("subset_mp", "full_mp")) |> 
  dplyr::left_join(jaccard_df, by = c("subset_mp", "full_mp")) |> 
  dplyr::mutate(
    subset_id = opt$subset_id,
    optimal_k = subset_mp_obj$n_metaprograms
  ) |> 
  dplyr::relocate(
    subset_id,
    optimal_k,
    .before = subset_mp
  )

readr::write_tsv(all_stats_df, opt$output_file)
