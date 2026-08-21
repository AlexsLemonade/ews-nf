#!/usr/bin/env Rscript

# This script downloads the Hallmark and Reactome gene sets from MSigDB using `msigdbr`
# and saves them as a single TSV file, `genesets/hallmark-reactome-gene-sets.tsv`.
# This script is not part of the workflow and only needs to be re-run if the gene sets need to be updated.

output_file <- here::here("genesets", "hallmark-reactome-gene-sets.tsv")

# Hallmark gene sets
hallmark_df <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
) |>
  dplyr::mutate(collection = "Hallmark") |>
  dplyr::select(collection, gs_name, ensembl_gene)

# Reactome gene sets
reactome_df <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:REACTOME"
) |>
  dplyr::mutate(collection = "Reactome") |>
  dplyr::select(collection, gs_name, ensembl_gene)

# combine and export
hallmark_reactome_df <- dplyr::bind_rows(hallmark_df, reactome_df) |>
  unique()

readr::write_tsv(hallmark_reactome_df, output_file)
