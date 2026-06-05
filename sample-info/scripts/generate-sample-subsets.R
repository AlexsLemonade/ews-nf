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

#' Use the number of existing reps and the total desired reps to figure 
#' out how many reps to create
#'
#' @param size Number of samples to include in the replicate
#' @param existing_counts List of existing sample counts for each subset
#'
#' @returns Number of replicates needed 
get_num_replicates <- function(size, existing_counts){
  # how many subsets of this size already exist
  n_existing <- existing_counts |>
    dplyr::filter(size == !!size) |>
    dplyr::pull(n)
  
  if(length(n_existing) == 0){
    n_existing <- 0
  }
  
  n_needed <- opt$n_replicates - n_existing 
  
  # print out a message about skipping the subsets if we already have enough 
  if (n_needed <= 0) {
    message(glue::glue("Skipping size {size}: already have {n_existing} subsets"))
  }
  
  return(n_needed)
}



# set up -----------------------------------------------------------------------

set.seed(2026)

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--n_replicates"),
    type = "numeric",
    default = 10,
    help = "Number of replicates for each subset size"
  ),
  make_option(
    opt_str = c("--max_patient_fraction"),
    type = "numeric",
    default = 0.1,
    help = "Maximum fraction of samples that can be from a single patient in a single subset"
  ),
  make_option(
    opt_str = c("--existing_subsets_file"),
    type = "character",
    default = here::here("sample-info", "sample-subsets.tsv"),
    help = "Path to existing subsets TSV file to append to"
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
subset_sizes <- seq(10, 90, 10)

# Check for existing subsets ---------------------------------------------------

# read in existing subsets if provided
if (!is.null(opt$existing_subsets_file)) {
  existing_subset_df <- readr::read_tsv(opt$existing_subsets_file, 
                                        col_types = list(.default = "c"))
  
  # get the last subset number from the full file BEFORE filtering
  # we don't want to overwrite existing IDs which might get confusing with processing
  last_subset_n <- existing_subset_df |>
    dplyr::pull(subset_id) |>
    stringr::str_extract("[0-9]+$") |>
    as.integer() |>
    max()
  
  # now filter to only the sizes we care about
  existing_subset_df <- existing_subset_df |>
    dplyr::mutate(subset_size = dplyr::n_distinct(scpca_sample_id), .by = subset_id) |>
    dplyr::filter(subset_size %in% subset_sizes)
  
  # extract the sets of sample ids already used, for duplicate checking
  existing_subsets_list <- existing_subset_df |>
    dplyr::summarise(samples = list(sort(scpca_sample_id)), .by = "subset_id") |>
    dplyr::mutate(subset_size = purrr::map_int(samples, length)) |>
    dplyr::summarise(subsets = list(samples), .by = "subset_size") |>
    tibble::deframe()
  
  # count how many subsets already exist per size
  existing_counts <- existing_subset_df |>
    dplyr::mutate(size = dplyr::n_distinct(scpca_sample_id), .by = subset_id) |>
    dplyr::distinct(subset_id, size) |>
    dplyr::count(size)
} else {
  # otherwise make these empty 
  existing_subsets_list <- list()
  last_subset_n <- 0
  existing_counts <- data.frame(size = integer(), n = integer())
}

# Get the subsets --------------------------------------------------------------

# data frame of all possible subsets 
subset_df <- subset_sizes |> 
  purrr::map(\(size){
    
    # get the number of reps needed 
    n_needed <- get_num_replicates(size, existing_counts)
    
    if(n_needed > 0) {
      
      # set starting point where duplicates are present
      # this ensures the first subset is grabbed and checked for duplicates 
      duplicates_present <- TRUE
      
      while(any(duplicates_present)){
        # get a dataframe with one column per replicate of subsets 
        subsets <- replicate(
          n_needed, 
          get_sample_subset(sample_df, size, opt$max_patient_fraction)
        )
        
        # pull out new subsets and existing subsets for duplicate checks
        new_subset_lists <- asplit(subsets, 2) |> 
          purrr::map(sort)
        
        existing_list <- existing_subsets_list[[as.character(size)]] |>
          purrr::map(sort)
        
        # check that there's no identical subsets that were chosen
        # if so, try again
        duplicates_present <- c(
          new_subset_lists,
          existing_list
        ) |> 
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
          subset_id = glue::glue("{size}-{replicate}"),
          subset_size = size
        )
    } else {
      return(NULL)
    }
  }) |> 
  dplyr::bind_rows()

full_cohort_subsets <- data.frame(
  scpca_sample_id = rep(sample_df$scpca_sample_id, 10),
  subset_id = rep(paste0("full-V", 1:10), each = nrow(sample_df)),
  subset_size = length(sample_df$scpca_sample_id)
)

# Export -----------------------------------------------------------------------

new_subset_df <- dplyr::bind_rows(subset_df, full_cohort_subsets) |> 
  dplyr::mutate(
    # reassign subset ids as subset_00X
    # start from the previously existing subset id
    subset_id = sprintf("subset_%03d", last_subset_n + as.integer(factor(subset_id, levels = unique(subset_id))))
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
    participant_id,
    subset_size
  ) 


# bind with existing and export
if (!is.null(opt$existing_subsets_file)) {
  dplyr::bind_rows(existing_subset_df, new_subset_df) |>
    readr::write_tsv(output_file)
} else {
  new_subset_df |>
    readr::write_tsv(output_file)
}
