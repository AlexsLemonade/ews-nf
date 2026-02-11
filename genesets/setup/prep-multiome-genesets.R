#!/usr/bin/env Rscript

# This script is used to prepare the multiome-module-genes.tsv file 
# the input file Furlan_multiome_3_subgroup_links.csv was originally generated from 
# https://github.com/furlan-lab/EwS_multiome/blob/669d6ef60a7946ece289ecb5ba353e7fa2c8a2c3/code/fig2.Rmd#L115-L136
# The output includes a gene_symbol, ensembl_gene_id, kmeans_group, and cell_state_module for each gene 
# only unique genes for each module are included in the output

# The gene sets are then added as new rows to the existing custom-gene-sets.tsv


multiome_original_file <- here::here("genesets", "Furlan_multiome_3_subgroup_links.csv")
multiome_output_file <- here::here("genesets", "multiome-module-genes.tsv")

# custom gene set file to add multiome list to 
custom_geneset_file <- here::here("genesets", "custom-gene-sets.tsv")

# Prep multiome gene sets ------------------------------------------------------

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

# get ensembl 110 which matches the version used for cellranger
ensembl_110 <- biomaRt::useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  version = 110
)

# map gene sybmols to ensembl gene id
gene_map_df <- biomaRt::getBM(
  attributes =  c("hgnc_symbol", "ensembl_gene_id"), 
  filters = "hgnc_symbol", 
  values = multiome_genesets_df$gene_symbol, 
  mart = ensembl_110
)

# join in ensembl gene ids
multiome_genesets_df <- multiome_genesets_df |> 
  dplyr::left_join(gene_map_df, by = c("gene_symbol" = "hgnc_symbol"))

# export new file with gene ids
readr::write_tsv(multiome_genesets_df, multiome_output_file)

# add multiome gene sets to custom gene set file -------------------------------
# read in existing custom gene sets
custom_genesets_df <- readr::read_tsv(custom_geneset_file)

# select only the columns we want to add 
# rename multiome_module as geneset
multiome_genesets_df <- multiome_genesets_df |> 
  dplyr::select(
    "geneset" = multiome_module,
    gene_symbol, 
    ensembl_gene_id
  )

# combine with existing custom list
all_genesets_df <- dplyr::bind_rows(
  custom_genesets_df,
  multiome_genesets_df
) |> 
  # in case we ever run this again, we don't want to duplicate the entries
  unique()

readr::write_tsv(all_genesets_df, custom_geneset_file)
