#!/usr/bin/env Rscript

# this script is used to run over representation analysis on a given set of metaprograms 
# the default output produces the manuscript/full-cohort-metaprograms/k-7_ora-results.tsv file

# functions and code below are copied from 04-metaprogram-metrics.R
# https://github.com/AlexsLemonade/ews-nf/blob/main/modules/metaprograms/resources/usr/bin/04-metaprogram-metrics.R

# function for generating a term2gene df with specific MSigDB gene sets
# input dataframe should include the geneset name, category and subcategory
# output includes a df with gs_name and ensembl_gene
prep_msigdb_genesets <- function(msigdb_df){
  
  # get list of categories that we need to grab from msigdb
  msigdb_category_df <- msigdb_df |>
    dplyr::select(category, subcategory) |>
    unique() 
  
  # list of msig genesets names
  # we have to pull all genesets from a specific category/subcategory 
  # and then filter to only the ones we care about 
  geneset_list <- msigdb_df$geneset |>
    purrr::set_names(msigdb_df$name)
  
  # first pull out info for each category and then pull out specific genes for geneset
  msigdb_geneset_df <- msigdb_category_df |>
    purrr::pmap(\(category, subcategory){
      
      # ifelse complains about null so need to use a full if statement 
      # to replace missing subcategory with NULL 
      if(is.na(subcategory)){
        subcategory <- NULL
      } else {
        subcategory <- subcategory
      }
      
      msigdbr::msigdbr(
        species = "Homo sapiens",
        collection = category,
        subcollection = subcategory
      )
    }) |>
    dplyr::bind_rows() |>
    # only keep relevant gene sets
    dplyr::filter(gs_name %in% geneset_list) |> 
    dplyr::select(gs_name, ensembl_gene)
  
  return(msigdb_geneset_df)
  
}


# function to extract the top genes from each MP 
# returns a list of vectors with the genes for each MP
# grabs top 200 genes by default
extract_top_genes <- function(mp_list, num_genes = 200){
  
  # get the top X number of genes for each MP
  mp_top_list <- mp_list |> 
    purrr::map(\(weights){
      weights |> 
        sort(decreasing = TRUE) |> 
        head(n=num_genes) |> 
        names() # return gene ids not weights
    })
  
  return(mp_top_list)
}

# function for running ORA with clusterProfiler::enricher on each metaprogram
# uses the top 200 genes for each metaprogram by default 
# returns a df with all gene sets for each metaprogram with p.adjust <= 0.05
run_ora <- function(mp_top_list, gene_universe, term2gene_df){
  
  # run ORA using the top 200 genes
  ora_results <- mp_top_list |> 
    purrr::map(\(genes){
      clusterProfiler::enricher(
        gene = genes,
        universe = gene_universe,
        TERM2GENE = term2gene_df
      )
    })
  
  # combine ora results into a single dataframe
  # only include genesets with pvalue < 0.05, the default by enricher 
  ora_results_df <- ora_results |> 
    purrr::map(\(ora) ora@result ) |> 
    purrr::list_rbind(names_to = "metaprogram") |>
    dplyr::filter(p.adjust <= 0.05) 
  
  return(ora_results_df)
  
}

# Read in arguments ------------------------------------------------------------

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--metaprograms_file"),
    type = "character",
    default = here::here("ews-nf-out/metaprograms/metaprograms-7_cnmf-5-10-15_filtered/metaprograms.rds"),
    help = "Path to RDS file containing metaprograms"
  ),
  make_option(
    opt_str = c("--custom_geneset_file"),
    type = "character",
    default = here::here("genesets/custom-gene-sets.tsv"),
    help = "Path to TSV file with custom gene sets to use with ORA"
  ),
  make_option(
    opt_str = c("--msigdb_geneset_file"),
    type = "character",
    default = here::here("genesets/msigdb-gene-sets.tsv"),
    help = "Path to TSV file with MSigDB gene sets to use with ORA"
  ),
  make_option(
    opt_str = c("--num_genes"),
    type = "integer",
    default = 200,
    help = "Number of genes to include with ORA. Genes are ranked based on gene weight and the top X genes are used. 
      Default: 200"
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    default = here::here("manuscript/full-cohort-metaprograms/k-7_ora-results.tsv"),
    help = "Path to TSV file to save the ORA results"
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = 2026,
    help = "A random seed for reproducibility."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Main script ------------------------------------------------------------------

# set the seed
set.seed(opt$seed)

# make sure inputs are correct
stopifnot(
  "Metaprograms file does not exist" = file.exists(opt$metaprograms_file),
  "Custom gene sets file does not exist" = file.exists(opt$custom_geneset_file),
  "MSigDB gene sets file does not exist" = file.exists(opt$msigdb_geneset_file),
  "Output file must end in .rds" = stringr::str_detect(opt$output_file, "\\.tsv$")
)

# object with metaprogram results
mp_object <- readr::read_rds(opt$metaprograms_file)

# read in custom gene lists and EF msigdb gene sets to use for ORA
ef_msigdb_df <- readr::read_tsv(opt$msigdb_geneset_file)
custom_geneset_df <- readr::read_tsv(opt$custom_geneset_file) |> 
  # reformat column names to match how msig pulls them
  dplyr::select(gs_name = geneset, ensembl_gene = ensembl_gene_id)

# grab hallmarks and reactome to use in addition to the custom and Ews specific gene lists 
hallmarks_msigdb_df <- msigdbr::msigdbr(
  species = "Homo sapiens", 
  collection = "H"
) |> 
  dplyr::select(gs_name, ensembl_gene)
reactome_msigdb_df <- msigdbr::msigdbr(
  species = "Homo sapiens", 
  collection = "C2",
  subcollection = "CP:REACTOME"
) |> 
  dplyr::select(gs_name, ensembl_gene)

# get genes in the EwS specific gene sets from MSigDB
ef_msigdb_geneset_df <- prep_msigdb_genesets(ef_msigdb_df)

# create a single term2genedf for running ORA with custom, Ews MSig, and hallmarks
all_term_df <- ef_msigdb_geneset_df |>
  dplyr::bind_rows(custom_geneset_df) |> 
  dplyr::bind_rows(hallmarks_msigdb_df) |> 
  dplyr::bind_rows(reactome_msigdb_df) |> 
  tidyr::drop_na()

# extract top genes
top_genes <- extract_top_genes(mp_object$metaprogram_list, opt$num_genes)

# get ora results
ora_results_df <- run_ora(
  mp_top_list = top_genes,
  gene_universe = mp_object$gene_universe,
  term2gene_df = all_term_df
)

# export 
readr::write_tsv(ora_results_df, opt$output_file)
