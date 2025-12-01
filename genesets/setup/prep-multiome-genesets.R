#!/usr/bin/env Rscript

# This script is used to prepare the multiome-module-genes.tsv file 
# the input file Furlan_multiome_3_subgroup_links.csv was originally generated from 
# https://github.com/furlan-lab/EwS_multiome/blob/669d6ef60a7946ece289ecb5ba353e7fa2c8a2c3/code/fig2.Rmd#L115-L136
# The output includes a gene_sybmol, ensembl_gene_id, kmeans_group, and cell_state_module for each gene 
# only unique genes for each module are included in the output

library(AnnotationDbi)
library(org.Hs.eg.db)


multiome_original_file <- here::here("genesets", "Furlan_multiome_3_subgroup_links.csv")
multiome_output_file <- here::here("genesets", "multiome-module-genes.tsv")

# select just the module and gene symbol columns, rename for clarity 
# add in a useful name other than just the kmean value for downstream use and plots
multiome_genesets_df <- readr::read_csv(multiome_original_file) |> 
  dplyr::select(gene_symbol = gene, kmeans_group = kmean) |> 
  dplyr::mutate(
    multiome_module = glue::glue("ews_multiome_module_{kmeans_group}")
  )

# filter out any duplicated gene symbols across programs before converting and saving 
# original paper only uses unique genes in plots and only top 300 genes for ORA to annotate the modules 
multiome_genesets_df <- multiome_genesets_df[!duplicated(multiome_genesets_df$gene_symbol),]

# get gene symbols and convert to ensembl IDs
symbols <- multiome_genesets_df$gene_symbol
gene_ids <- mapIds(
  org.Hs.eg.db, 
  keys = symbols, 
  column = "ENSEMBL", 
  keytype = "SYMBOL", 
  multiVals = "first"
)

# check that everything is equal before adding it back 
stopifnot("Error converting all gene symbols" = all.equal(multiome_genesets_df$gene_symbol, names(gene_ids)))

# add new column and export 
multiome_genesets_df$ensembl_gene_id <- gene_ids
readr::write_tsv(multiome_genesets_df, multiome_output_file)
