#!/usr/bin/env Rscript

# This script is used to generate a GeneSetCollection containing a 
# group of custom and MSigDB genesets to be used for running AUCell on a processed SCE object
# The GeneSetCollection is saved as an RDS file

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--custom_genesets_file"),
    type = "character",
    help = "Path to TSV file containing a list of custom gene sets.
      Must contain a `geneset` column with the name of the gene set and `ensembl_gene_id` with the Ensembl gene identifier."
  ),
  make_option(
    opt_str = c("--msigdb_genesets"),
    type = "character",
    help = "Path to TSV file containing all gene sets from MSigDB to use with AUCell.
      Must contain columns with `name`, `geneset`, `category`, and `subcategory`."
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    help = "Path to RDS file where GeneSetCollection will be saved"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# make sure input files exist
stopifnot(
  "MSigDB gene set file does not exist" = file.exists(opt$msigdb_genesets),
  "Custom gene set file does not exist" = file.exists(opt$custom_genesets_file),
  "Output file must end in .rds" = stringr::str_detect(opt$output_file, "\\.rds$")
)

# make sure directory exists for writing output
output_dir <- dirname(opt$output_file)
fs::dir_create(output_dir)

# read in gene sets
msig_genesets_df <- readr::read_tsv(opt$msigdb_genesets)
custom_genesets_df <- readr::read_tsv(opt$custom_genesets_file)

# Msigdb genee sets ------------------------------------------------------------

# get list of categories that we need to grab from msigdb
category_list <- msig_genesets_df |>
  dplyr::select(category, subcategory) |>
  unique() |>
  purrr::transpose()

# list of genesets and names
geneset_list <- msig_genesets_df$geneset |>
  purrr::set_names(msig_genesets_df$name)

# pull gene sets from msigbdr
# first pull out info for each category and then pull out specific genes for geneset
msig_genes_df <- category_list |>
  purrr::map(\(category_list){
    # replace subcategory with default NULL
    # can't use NULL in tsv since it gets read in as a character
    if (is.na(category_list$subcategory)) {
      subcategory <- NULL
    } else {
      subcategory <- category_list$subcategory
    }
    
    msigdbr::msigdbr(
      species = "Homo sapiens",
      collection = category_list$category,
      subcollection = subcategory
    )
  }) |>
  dplyr::bind_rows() |>
  # only keep relevant gene sets
  dplyr::filter(gs_name %in% geneset_list)

# create named list of genes in each gene set 
msig_genes_list <- geneset_list |>
  purrr::map(\(name){
    genes <- msig_genes_df |>
      dplyr::filter(gs_name == name) |>
      dplyr::pull(ensembl_gene) |>
      unique()
  })

# Custom gene sets -------------------------------------------------------------

# get names of gene sets using name of the files
custom_geneset_names <- unique(custom_genesets_df$geneset)

# read in custom gene sets
custom_genes_list <- custom_geneset_names |>
  purrr::map(\(name) {
    gene_ids <- custom_genesets_df |> 
      dplyr::filter(geneset == name) |>
      dplyr::pull(ensembl_gene_id) |>
      unique()
  }) |> 
  purrr::set_names(custom_geneset_names)

# Build GeneSetCollection ------------------------------------------------------

# combine custom and msig
genes_list <- c(msig_genes_list, custom_genes_list)

# build GeneSetCollection for AUCell
collection <- genes_list |>
  purrr::imap(\(genes, name) GSEABase::GeneSet(genes, setName = name)) |>
  GSEABase::GeneSetCollection()

readr::write_rds(collection, opt$output_file)
