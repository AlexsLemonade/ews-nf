#!/usr/bin/env Rscript

# This script reads in the sample metadata file and creates random subsets of samples 
# of various sizes to use for identifying the point of saturation for metaprograms 

# The output is a TSV file with subset id, sample id and patient id

# Functions ----

#' Obtain a subset of samples from the full sample data frame
#'
#' @param sample_df Dataframe with scpca_sample_id and participant_id
#' @param size Size of subset to obtain, e.g. 10
#'
#' @returns A list of sample ids with length equivalent to `size`
get_sample_subset <- function(sample_df, size, max_patient_fraction){
  
  # set the per patient cap as a fraction of all samples
  per_patient_cap <- floor(max_patient_fraction * size)
  # start with this value being above the per patient cap 
  max_samples_per_patient <- per_patient_cap + 1
  
  # as long as the max samples per patient > cap get a subset 
  while(max_samples_per_patient > per_patient_cap) {
    
    sample_subset <- sample(sample_df$scpca_sample_id, size)
    
    max_samples_per_patient <- sample_df |> 
      dplyr::filter(scpca_sample_id %in% sample_subset) |> 
      dplyr::count(participant_id) |> 
      dplyr::pull(n) |> 
      max() 
  }
  
  # sort so we can check that we don't get any duplicate subsets
  return(sort(sample_subset))
  
}

# set up -----------------------------------------------------------------------

set.seed(2026)

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--n_replicates"),
    type = "numeric",
    default = 5,
    help = "Number of replicates for each subset size"
  ),
  make_option(
    opt_str = c("--max_patient_fraction"),
    type = "numeric",
    default = 0.1,
    help = "Maximum fraction of samples that can be from a single patient in a single subset"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))


# define intput and output files
sample_metadata_file <- here::here("sample-info", "sample-metadata.tsv")
blacklist_sample_file <- here::here("sample-info", "low-cell-count-libraries.txt")

output_file <- here::here("sample-info", "sample-subsets.tsv")

# read in the list of samples that shouldn't be included
# these samples have < 200 cells
blacklist_sample_ids <- readLines(blacklist_sample_file)

# read in sample dataframe and remove blacklisted samples
sample_df <- readr::read_tsv(sample_metadata_file, col_types = list(.default = "c")) |> 
  dplyr::select(scpca_sample_id, participant_id) |> 
  dplyr::filter(!scpca_sample_id %in% blacklist_sample_ids)

# define all the subset sizes we'll create 
subset_sizes <- seq(10, 100, 10)

# Get the subsets --------------------------------------------------------------

# data frame of all possible subsets 
subset_df <- subset_sizes |> 
  purrr::map(\(size){
    
    # set starting point where duplicates are present
    # this ensures the first subset is grabbed and checked for duplicates 
    duplicates_present <- TRUE
    
    while(any(duplicates_present)){
      # get a dataframe with one column per replicate of subsets 
      subsets <- replicate(
        opt$n_replicates, 
        get_sample_subset(sample_df, size, opt$max_patient_fraction)
      )
      
      # check that there's no identical subsets that were chosen
      # if so, try again
      duplicates_present <- t(subsets) |> 
        duplicated()
    }
    
    # convert to one row per sample id per subset
    subsets |> 
      as.data.frame() |> 
      tidyr::pivot_longer(cols = everything(), 
                          names_to = "replicate",
                          values_to = "scpca_sample_id") |> 
      dplyr::mutate(
        # assign a temporary subset id
        subset_id = glue::glue("{size}-{replicate}")
      )
  }) |> 
  dplyr::bind_rows()

# Export -----------------------------------------------------------------------

subset_df |> 
  dplyr::mutate(
    # reassign subset ids as subset_00X
    subset_id = sprintf("subset_%03d", as.integer(factor(subset_id, levels = unique(subset_id))))
  ) |> 
  # join in participant id information
  dplyr::left_join(
    sample_df,
    by = "scpca_sample_id"
  ) |> 
  # reorder and select specific columns
  dplyr::select(
    subset_id,
    scpca_sample_id, 
    participant_id
  ) |> 
  readr::write_tsv(output_file)

