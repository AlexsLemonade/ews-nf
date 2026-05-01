#!/usr/bin/env Rscript

# This script creates the list of reference cell types used by the workflow to 
# use as a normal pooled reference with AUCell and InferCNV 

# input files 
# list of validation groups to keep in the normal reference
ref_validation_groups_file <- here::here("references", "reference-validation-groups.txt")
# all possible cell types and validation groups
consensus_validation_groups_file <- here::here("references", "consensus-validation-groups.tsv")

# output file
ref_celltypes_file <- here::here("references", "reference-celltypes.tsv")

# read in validation groups to allow
ref_validation_groups <- readLines(ref_validation_groups_file)

# read in all cell types and filter 
ref_celltypes_df <- readr::read_tsv(consensus_validation_groups_file) |> 
  dplyr::filter(validation_group_annotation %in% ref_validation_groups) |> 
  dplyr::select(consensus_annotation, consensus_ontology)

readr::write_tsv(ref_celltypes_df, ref_celltypes_file)


old_ref_celltypes <- readr::read_tsv("references/reference-celltypes.tsv") |> 
  dplyr::pull(consensus_annotation)
