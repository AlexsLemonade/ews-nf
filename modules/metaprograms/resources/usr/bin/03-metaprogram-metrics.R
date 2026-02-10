#!/usr/bin/env Rscript

# This script is used to calculate a set of metrics on metaprograms
# The output is a RDS file with 

# mp_metrics_df: A dataframe of metrics calculated for each MP
# geneset_prop_df: A dataframe with the proportion of total MPs that each 
#  gene set is associated with as determined by ORA
# ora_results_df: A dataframe of all results output from ORA

# mp_metrics_df contains the following metrics: 

# max_correlation: Max correlation of the specified metaprogram to all other metaprograms
# median_difference: Difference in medians between the correlation of spectra 
#  within metaprograms and between metaprograms
# num_samples: Number of samples represented in the metaprogram
# proportion_of_samples: Proportion of the total samples represented in the metaprogram
# num_sig_genesets: Total number of gene sets identified as significant by ORA
# max_jaccard: Max Jaccard index for each metaprogram calculated by comparing 
#  the list of gene sets for each metaprogram after ORA
#  this column includes all gene sets
# max_jaccard_no_ews: Max Jaccard index after removing EF-activated and EF-repressed gene sets

# Additional entries in the object will include information about the parameters
# used for metaprogram generation 

# Functions --------------------------------------------------------------------

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
# returns a nested list with the genes for each MP
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
run_ora <- function(mp_top_list, gene_universe, term2gene_df, num_genes=200){
  
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

# helper function for getting jaccard between two vectors
calc_jaccard <- function(vec1, vec2) {
  length(intersect(vec1, vec2)) / length(union(vec1, vec2))
}

# get the max jaccard for overlap between two groups for a set of metaprograms
# groups_to_compare is a nested list of items to compare with names of the list corresponding to MPs 
# returns a df with the off diagonal max jaccard for each metaprogram
max_jaccard <- function(groups_to_compare){

  # create the grid of comparisons
  cross_df <- tidyr::expand_grid(id1 = names(groups_to_compare), id2 = names(groups_to_compare))
  
  jaccard_df <- cross_df |>
    # remove comparisons to self
    dplyr::filter(id1 != id2) |> 
    dplyr::rowwise() |> 
    # calculate a single Jaccard index for each combination using split lists & ids
    dplyr::mutate(jaccard = calc_jaccard(groups_to_compare[[id1]], groups_to_compare[[id2]])) |>
    dplyr::ungroup() |> # remove rowwise prior to getting max
    # get the max jaccard
    dplyr::slice_max(jaccard, by = id1, with_ties = FALSE) |>
    dplyr::rename(metaprogram = id1, max_jaccard = jaccard) |> 
    dplyr::select(metaprogram, max_jaccard)
  
  return(jaccard_df)
}

# Read in arguments ------------------------------------------------------------

library(optparse)
library(SingleCellExperiment)

option_list <- list(
  make_option(
    opt_str = c("--metaprograms_file"),
    type = "character",
    default = "metaprograms.rds",
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
    opt_str = c("--category_file"),
    type = "character",
    default = here::here("genesets/gene-set-categories.tsv"),
    help = "Path to TSV file with categories for all gene sets"
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    default = "mp_metrics.rds", 
    help = "Path to RDS file to save the metaprogram and gene set metrics"
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
  "Gene set categories file does not exist" = file.exists(opt$category_file),
  "Output file must end in .rds" = stringr::str_detect(opt$output_file, "\\.rds$")
)

# object with metaprogram results
mp_object <- readr::read_rds(opt$metaprograms_file)

# read in custom gene lists and EF msigdb gene sets to use for ORA
ef_msigdb_df <- readr::read_tsv(opt$msigdb_geneset_file)
ef_category_df <- readr::read_tsv(opt$category_file)
custom_geneset_df <- readr::read_tsv(opt$custom_geneset_file) |> 
  # reformat column names to match how msig pulls them
  dplyr::select(gs_name = geneset, ensembl_gene = ensembl_gene_id)

# Max correlation --------------------------------------------------------------

# combine all mps into a matrix and get correlation to each other
mp_mtx <- do.call(cbind, mp_object$metaprogram_list)
mp_cor_matrix <- cor(mp_mtx, method = "pearson", use = "pairwise.complete.obs")

# the diagonals are all self comparisons so set those to 0 
diag(mp_cor_matrix) <- 0
max_cor_df <- data.frame(
  metaprogram = rownames(mp_cor_matrix), 
  max_correlation = rowMax(mp_cor_matrix)
)

# Difference in medians --------------------------------------------------------

# convert the mp cluster list into a dataframe with one row per spectra
# use this to get the assigned mp for each cluster/metaprogram when calculating metrics
mp_lookup <- tibble::enframe(mp_object$cluster_list, name = "metaprogram", value = "spectra") |> 
  tidyr::unnest_longer(spectra)

# for each correlation, add in the MPs being compared and classify MP
mp_cor_df <- mp_object$spectra_cor_df |> 
  dplyr::filter(spectra1 != spectra2) |> 
  # add in the mp assignment for the first spectra
  dplyr::left_join(mp_lookup, by = c("spectra1" = "spectra")) |> 
  dplyr::rename(mp1 = metaprogram) |> 
  # add in the metaprogram assignment for the second spectra
  dplyr::left_join(mp_lookup, by = c("spectra2" = "spectra")) |> 
  dplyr::rename(mp2 = metaprogram) |>
  # add a column indicating if the spectra are assigned to the same or different mps
  dplyr::mutate(within_mp = ifelse(mp1 == mp2, "within", "between"))

# calculate the difference in medians between the two groups to display on the plot
mp_medians_df <- mp_cor_df |> 
  dplyr::group_by(mp1, within_mp) |> 
  dplyr::summarise(
    group_median = median(correlation)
  ) |> 
  tidyr::pivot_wider(
    id_cols = mp1,
    names_from = within_mp,
    values_from = group_median
  ) |> 
  dplyr::mutate(
    median_difference = within - between
  ) |> 
  dplyr::select(metaprogram = mp1, median_difference)

# Sample proportion ------------------------------------------------------------

# get the total number of samples to use as the denominator for proportion
# use all samples even if they are not all present in the spectra after filtering
total_samples <- length(mp_object$unique_ids)

# count the total number and proportion of samples per mp
mp_prop_df <- mp_lookup |>
  dplyr::mutate(sample_id = stringr::str_extract(spectra, "SCPCS[0-9]{6}")) |>
  dplyr::select(metaprogram, sample_id) |> 
  dplyr::distinct() |> 
  dplyr::group_by(metaprogram) |>
  dplyr::summarize(
    num_samples = length(sample_id),
    proportion_of_samples = num_samples / total_samples
  )

# Overlap in top genes ---------------------------------------------------------

# get list of top genes for each mp 
top_genes <- extract_top_genes(mp_list = mp_object$metaprogram_list)

# get the max jaccard between top 200 genes for each MP 
mp_genes_overlap_df <- max_jaccard(top_genes) |> 
  dplyr::rename(top_genes_max_jaccard = max_jaccard)


# ORA --------------------------------------------------------------------------

# grab hallmarks to use in addition to the custom and Ews specific gene lists 
hallmarks_msigdb_df <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") |> 
  dplyr::select(gs_name, ensembl_gene)

# get genes in the EwS specific gene sets from MSigDB
ef_msigdb_geneset_df <- prep_msigdb_genesets(ef_msigdb_df)

# create a single term2genedf for running ORA with custom, Ews MSig, and hallmarks
all_term_df <- ef_msigdb_geneset_df |>
  dplyr::bind_rows(custom_geneset_df) |> 
  dplyr::bind_rows(hallmarks_msigdb_df) |> 
  tidyr::drop_na()

# get ora results
ora_results_df <- run_ora(
  mp_top_list = top_genes,
  gene_universe = mp_object$gene_universe,
  term2gene_df = all_term_df
)

# table of number of gene sets per metaprogram 
num_genesets_df <- ora_results_df |> 
  dplyr::group_by(metaprogram) |> 
  dplyr::summarise(num_sig_genesets = length(ID))

# Gene set proportion ---------------------------------------------------------_

# table of genesets and the proportion of MPs 
geneset_prop_df <- ora_results_df |> 
  dplyr::group_by(ID) |>
  dplyr::summarise(
    # get the num mps and proportion that each gene set is significant in
    num_mps = dplyr::n_distinct(metaprogram), 
    proportion_mps = num_mps/mp_object$n_metaprograms,
    metaprograms = list(metaprogram)
  ) |> 
  dplyr::arrange(desc(proportion_mps)) |> 
  as.data.frame()

# Jaccard index ----------------------------------------------------------------

# max jaccard between all gene sets identified for each MP
# get a nested list of groups to compare with names as the MPs and gene sets for each MP
geneset_groups <- split(ora_results_df$ID, ora_results_df$metaprogram)
max_jaccard_df <- max_jaccard(groups_to_compare = geneset_groups)

# get the EF activated and EF repressed gene sets to remove from the results
ews_fli1_genesets <- ef_category_df |> 
  dplyr::filter(geneset_category %in% c("EF-activated", "EF-repressed")) |> 
  dplyr::pull(full_geneset_name)

# filter out EF gene sets and recalculate jaccard
filtered_ora_results_df <- ora_results_df |> 
  dplyr::filter(!(ID %in% ews_fli1_genesets))
filtered_geneset_groups <- split(filtered_ora_results_df$ID, filtered_ora_results_df$metaprogram)

max_jaccard_no_ews <- max_jaccard(filtered_geneset_groups) |> 
  dplyr::rename(max_jaccard_no_ews = max_jaccard)

# Combine and export -----------------------------------------------------------

# join all mp metrics 
mp_metrics_df <- max_cor_df |> 
  dplyr::left_join(mp_medians_df, by = "metaprogram") |> 
  dplyr::left_join(mp_prop_df, by = "metaprogram") |> 
  dplyr::left_join(mp_genes_overlap_df, by = "metaprogram") |> 
  dplyr::left_join(num_genesets_df, by = "metaprogram") |>
  dplyr::left_join(max_jaccard_df, by = "metaprogram") |> 
  dplyr::left_join(max_jaccard_no_ews, by = "metaprogram") 

# create an object for exporting with mp metrics, gene set entropy, and ora results
metrics_obj <- list(
  mp_metrics_df = mp_metrics_df,
  geneset_prop_df = geneset_prop_df,
  ora_results_df = ora_results_df,
  # add some key parameters which will help when we combine these TSV files across different parameters
  n_metaprograms = mp_object$n_metaprograms,
  cnmf_k_range = mp_object$cnmf_k_range,
  filter_spectra = mp_object$filter_spectra,
  orphan_cutoff = mp_object$orphan_cutoff,
  removed_unique_ids = mp_object$removed_unique_ids
)

readr::write_rds(metrics_obj, opt$output_file)

