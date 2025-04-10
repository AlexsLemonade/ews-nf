#!/usr/bin/env Rscript

# This script is used to convert the gene symbol column of custom-gene-sets.tsv to Ensembl Ids
# The ensembl Ids are then saved to a column, `ensembl_gene_id`

library(AnnotationDbi)
library(org.Hs.eg.db)

gene_list_file <- here::here("genesets", "custom-gene-sets.tsv")
gene_df <- readr::read_tsv(gene_list_file)

# get gene symbols and convert to ensembl IDs
symbols <- gene_df$gene_symbol
gene_ids <- mapIds(
  org.Hs.eg.db, 
  keys = symbols, 
  column = "ENSEMBL", 
  keytype = "SYMBOL", 
  multiVals = "first"
)

# check that everything is equal before adding it back 
stopifnot("Error converting all gene symbols" = all.equal(gene_df$gene_symbol, names(gene_ids)))

# add new column and export 
gene_df$ensembl_gene_id <- gene_ids
readr::write_tsv(gene_df, gene_list_file)
