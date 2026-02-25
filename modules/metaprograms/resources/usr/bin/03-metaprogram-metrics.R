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
# num_sig_genesets_pvalue: P value associated with permutation testing for number of observed gene sets after ORA
# num_sig_genesets_adj_pvalue: Adjusted pvalue associated with permutation testing for number of gene sets after ORA
# num_unique_genesets: Total number of gene sets uniquely mapped to a single metaprogram 
# proportion_unique_genesetss: Total proportion of gene sets uniquely mapped to a single metaprogram
# unique_genesets_pvolue: P value associated with permutation testing for unique gene sets
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

# helper function for getting jaccard between two vectors
calc_jaccard <- function(vec1, vec2) {
  length(intersect(vec1, vec2)) / length(union(vec1, vec2))
}

# get the max jaccard for overlap between two groups for a set of metaprograms
# e.g., a list of genes or a list of gene sets assigned to each metaprogram 
# groups_list is a list of vectors containing the groups assigned to each metaprogram
# returns a df with the off diagonal max jaccard for each metaprogram
max_jaccard <- function(groups_list){

  # create the grid of comparisons
  cross_df <- tidyr::expand_grid(id1 = names(groups_list), id2 = names(groups_list))
  
  jaccard_df <- cross_df |>
    # remove comparisons to self
    dplyr::filter(id1 != id2) |> 
    dplyr::rowwise() |> 
    # calculate a single Jaccard index for each combination using split lists & ids
    dplyr::mutate(jaccard = calc_jaccard(groups_list[[id1]], groups_list[[id2]])) |>
    dplyr::ungroup() |> # remove rowwise prior to getting max
    # get the max jaccard
    dplyr::slice_max(jaccard, by = id1, with_ties = FALSE) |>
    dplyr::rename(metaprogram = id1, max_jaccard = jaccard) |> 
    dplyr::select(metaprogram, max_jaccard)
  
  return(jaccard_df)
}

# helper function to get a table with the total number of genesets per metaprogram 
get_total_genesets <- function(mps, ora_results_df){
  
  # table of number of gene sets per metaprogram 
  num_genesets_df <- ora_results_df |> 
    dplyr::group_by(metaprogram) |> 
    dplyr::summarise(num_sig_genesets = length(ID))
  
  df <- data.frame(
    metaprogram = mps
  ) |> 
    dplyr::left_join(num_genesets_df, by = "metaprogram") |> 
    dplyr::mutate(
      num_sig_genesets = dplyr::if_else(is.na(num_sig_genesets), 0, num_sig_genesets)
    )
  
  return(df)
  
}


# function to randomly shuffle genes assigned to metaprograms, recalculate the top scores, 
# and then redo ORA and calculate stats
# mp_mtx is a matrix with rownames as genes and columns as metaprograms
# here we report both the num of significant gene sets and the max jaccard for each metaprogram
permute_ora <- function(mp_mtx, nreps, gene_universe, term2gene_df) {
  
  reps <- seq(1, nreps)
  
  background <- replicate(nreps, {
      
      # shuffle the genes (rownames) in the matrix 
      rownames(mp_mtx) <- sample(rownames(mp_mtx))
      
      # extract top 200 genes for each mp 
      # first turn matrix into a list of vector weights where the names are gene ids
      mp_list <- asplit(mp_mtx, MARGIN = 2)
      top_genes <- extract_top_genes(mp_list)
      
      # run ORA 
      ora_results_df <- run_ora(
        top_genes,
        gene_universe,
        term2gene_df
      )
      
      # table of number of gene sets per metaprogram 
      get_total_genesets(names(mp_list), ora_results_df) 
    }, simplify = FALSE) |> # return as a list not mtx
    dplyr::bind_rows(.id = "replicate")
  
  return(background)
  
}

# calculate significance for permutation testing
# reports a pvalue for number of gene sets returned from ORA
ora_permutation_significance <- function(observed_df, background_df, nreps){
  
  # combine the observed values with the background
  combined_df <- observed_df |> 
    dplyr::select(
      metaprogram,
      obs_num_sig_genesets = num_sig_genesets
    ) |> 
    dplyr::left_join(background_df, by = c("metaprogram")) |> 
    dplyr::select(-replicate)
  
  # calculate the pvalue for num genesets 
  pvalue_df <- combined_df |> 
    dplyr::mutate(
      # indicate which rows have background > obs
      greater_than_obs_num_genesets = num_sig_genesets > obs_num_sig_genesets
    ) |> 
    dplyr::group_by(metaprogram) |> 
    dplyr::summarize(
      # calculate p values and retain a column with the observed value
      num_sig_genesets = unique(obs_num_sig_genesets),
      num_sig_genesets_pvalue = (sum(greater_than_obs_num_genesets) +1)/(nreps + 1)
    ) |> 
    # add adjusted pvalue 
    dplyr::mutate(
      num_sig_genesets_adj_pvalue = p.adjust(num_sig_genesets_pvalue, method = "BH")
    )
  
  return(pvalue_df)
  
}

# helper function to calculate the total proportion and number of uniquely mapped gene sets
# geneset_list is a named list of gene sets mapped to each MP after ORA
prop_unique_mapped_genesets <- function(geneset_list){
  
  # count gene sets that appear only once across all metaprograms
  all_observed_genesets <- unlist(geneset_list) 
  total_unique_genesets <- length(unique(all_observed_genesets))
  
  num_unique_genesets <- sum(table(all_observed_genesets) == 1)
  proportion_unique_genesets <- num_unique_genesets/total_unique_genesets
  
  num_geneset_list <- list(
    num_unique_genesets = num_unique_genesets,
    proportion_unique_genesets= proportion_unique_genesets
  )
  
}

# shuffle gene sets assigned to each MP using ORA
# we shuffle geneset assignments ensuring no duplicates per mp
# and calculate total numer and proportion of unique gene sets
# one value per value of k
permute_unique_genesets <- function(ora_results_df, nreps) {
  
  # pool of all significant gene sets across all metaprograms
  all_genesets <- unique(ora_results_df$ID)
  total_genesets <- length(all_genesets)
  
  # get the size (number of sig gene sets) for each metaprogram
  mp_geneset_counts <- ora_results_df |>
    dplyr::group_by(metaprogram) |>
    dplyr::summarise(geneset_count = dplyr::n_distinct(ID), .groups = "drop")
  
  # get the unique geneset counts so we only compute one distribution per count
  unique_counts <- unique(mp_sizes$geneset_count)
  
  # mp names
  mp_names <- mp_sizes$metaprogram
  
  # for each unique size, generate a background distribution of max jaccard
  # this keeps things efficient because we reuse the same distribution for all mps
  # that have the same number of sig gene sets
  size_distributions <- purrr::map(unique_sizes, \(size) {
    
    replicate(nreps, {
      
      # sample without replacement from the full pool
      # ensuring no duplicate gene sets within a metaprogram
      sampled_genesets <- sample(all_genesets, size, replace = FALSE)
      
      data.frame(
        size = size,
        geneset = sampled_genesets
      )
    }, simplify = FALSE) |> # return as a dataframe not mtx
      dplyr::bind_rows(.id = "replicate")
    
  }) |> 
    purrr::set_names(unique_sizes) |>
    dplyr::bind_rows()
  
  # for each rep, calculate max jaccard across metaprograms using shuffled assignments
  background <- replicate(
    nreps,
    {
      # for this rep, reconstruct geneset groups per mp using their observed size
      # but with randomly sampled gene sets
      shuffled_groups <- mp_names |>
        purrr::map(\(mp) {
          mp_size <- mp_sizes |>
            dplyr::filter(metaprogram == mp) |>
            dplyr::pull(size)
          
          size_distributions |>
            dplyr::filter(size == mp_size) |>
            dplyr::slice_sample(n = 1) |>
            dplyr::pull(geneset)
        }) |>
        purrr::set_names(mp_names)
      
      # get unique gene sets as data frame
      unique_stats <- prop_unique_mapped_genesets(shuffled_groups) |> 
        as.data.frame()
      
    },
    simplify = FALSE
  ) |> 
    dplyr::bind_rows(.id = "replicate")
  
  
  return(background)
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
    opt_str = c("--nreps"),
    type = "integer",
    default = 1000, 
    help = "Number of replicates to perform for permutation testing"
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
num_genesets_df <- get_total_genesets(names(mp_object$metaprogram_list), ora_results_df)

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
max_jaccard_df <- max_jaccard(groups_list = geneset_groups)

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

# Permute ORA ------------------------------------------------------------------

# get the backround distribution for num genesets from ORA
ora_background_df <- permute_ora(
  mp_mtx, 
  nreps = opt$nreps, 
  gene_universe = mp_object$gene_universe, 
  term2gene_df = all_term_df
)

# calculate pvalues 
ora_pvalue_df <- ora_permutation_significance(num_genesets_df, ora_background_df, opt$nreps)

# Proportion of uniquely mapped gene sets --------------------------------------
observed_unique_genesets <- prop_unique_mapped_genesets(geneset_groups)

# background distribution for total unique gene sets
unique_geneset_background_df <- permute_unique_genesets(ora_results_df, opt$nreps) |> 
  dplyr::mutate(
    # indicate which rows have background > obs
    greater_than_obs = num_unique_genesets > observed_unique_genesets$num_unique_genesets
  )

# get a single pvalue 
unique_genesets_pvalue <- (sum(unique_geneset_background_df$greater_than_obs) +1)/(opt$nreps + 1)

# Combine and export -----------------------------------------------------------

# join all mp metrics 
mp_metrics_df <- max_cor_df |> 
  dplyr::left_join(mp_medians_df, by = "metaprogram") |> 
  dplyr::left_join(mp_prop_df, by = "metaprogram") |> 
  dplyr::left_join(mp_genes_overlap_df, by = "metaprogram") |> 
  # includes total num genesets and associated pvalues 
  dplyr::left_join(ora_pvalue_df, by = "metaprogram") |> 
  # now add in the total proportion of uniquely mapped gene sets and pvalue 
  # only one value per set of metaprograms so all rows have the same value
  dplyr::mutate(
    num_unique_genesets = observed_unique_genesets$num_unique_genesets,
    proportion_unique_genesets = observed_unique_genesets$proportion_unique_genesets,
    unique_genesets_pavlue = unique_genesets_pvalue
  ) |>
  dplyr::left_join(max_jaccard_no_ews, by = "metaprogram") 

# create an object for exporting with mp metrics, gene set entropy, and ora results
metrics_obj <- list(
  mp_metrics_df = mp_metrics_df,
  geneset_prop_df = geneset_prop_df,
  ora_results_df = ora_results_df,
  ora_background_df = ora_background_df, # add in background stats for plotting
  unique_geneset_background_df = unique_geneset_background_df,
  # add some key parameters which will help when we combine these TSV files across different parameters
  n_metaprograms = mp_object$n_metaprograms,
  cnmf_k_range = mp_object$cnmf_k_range,
  filter_spectra = mp_object$filter_spectra,
  orphan_cutoff = mp_object$orphan_cutoff,
  removed_unique_ids = mp_object$removed_unique_ids
)

readr::write_rds(metrics_obj, opt$output_file)

