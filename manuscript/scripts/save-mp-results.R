#!/usr/bin/env Rscript

# This script reads in the metaprograms and extracts the top 200 genes for each MP
# and the full ORA results and saves them to TSV files in the repo 

# Before running this script, you will need to have the following files 
# copied from pitt's cluster (/vast/acillo/alhawk01)
# ews-nf-out/metaprograms/metaprograms-7_cnmf-5-10-15_filtered/k-7_metaprograms-metrics.rds
# scpca-nf-out/SCPCP000026/SCPCS000841/SCPCL001228-SCPCS000841_processed.rds

# input files 
metrics_file <- here::here("ews-nf-out", "metaprograms", "metaprograms-7_cnmf-5-10-15_filtered", "k-7_metaprograms-metrics.rds")
# use sce to get gene symbols in rowdata
sce_file <- here::here("scpca-nf-out", "SCPCP000026", "SCPCS000841", "SCPCL001228-SCPCS000841_processed.rds")

# output files 
output_mp_genes_file <- here::here("manuscript", "full-cohort-metaprograms", "k-7_metaprogram-top-genes.tsv")

# read in files 
mp_metrics_obj <- readr::read_rds(metrics_file)
sce <- readr::read_rds(sce_file)

# pull out the ids and gene symbols from the sce
id_to_symbol_df <- rowData(sce) |> 
  as.data.frame() |> 
  dplyr::select(gene_ids, gene_symbol)

# turn the list of top 200 genes for each mp into a data frame
mp_top_df <- mp_metrics_obj$top_genes |> 
  stack() |> 
  dplyr::rename(
    "metaprogram" = ind,
    "gene_ids" = values 
  ) |>
  dplyr::left_join(id_to_symbol_df, by = "gene_ids") |>
  # now just reorder
  dplyr::select(metaprogram, gene_symbol, gene_ids)

# export   
readr::write_tsv(mp_top_df, output_mp_genes_file)
