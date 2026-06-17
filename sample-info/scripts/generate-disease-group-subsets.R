#!/usr/bin/env Rscript

# This script reads in the sample metadata file and creates subsets 
# of all primary and all metastatic samples

# The output is a TSV file with subset id and sample id 

# define intput and output files
sample_metadata_file <- here::here("sample-info", "sample-metadata.tsv")
output_file <- here::here("sample-info", "disease-relevant-subsets.tsv")

# read in sample dataframe and remove blacklisted samples
sample_df <- readr::read_tsv(sample_metadata_file, col_types = list(.default = "c")) |> 
  dplyr::select(scpca_sample_id, tumor_type)

# primary samples 
primary_df <- sample_df |> 
  dplyr::filter(tumor_type == "primary") |> 
  dplyr::mutate(
    subset_id = "primary"
  )

# metastatic samples 
mets_df <- sample_df |> 
  dplyr::filter(tumor_type == "met") |> 
  dplyr::mutate(
    subset_id = "met"
  )

# combine and export 
dplyr::bind_rows(list(
  primary_df,
  mets_df
)) |> 
  dplyr::relocate(subset_id, .before = 1) |> 
  dplyr::select(-tumor_type) |> 
  readr::write_tsv(output_file)
