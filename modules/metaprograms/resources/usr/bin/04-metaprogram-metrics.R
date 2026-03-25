#!/usr/bin/env Rscript

# This script is used to calculate a set of metrics on metaprograms
# The output is a RDS file with:
# mp_metrics_df: A dataframe of metrics calculated for each MP
# ora_results_df: A dataframe of all results output from ORA
# background_mp_stats_df: A dataframe with the background values for metaprogram stats
# background_specificity_df: Dataframe with background for gene set specificity

# mp_metrics_df contains the following metrics: 

# metaprogram specific metrics ------
# max_correlation: Max correlation of the specified metaprogram to all other metaprograms
# median_difference: Difference in medians between the correlation of spectra 
#  within metaprograms and between metaprograms
# num_spectra_per_mp: Number of spectra found in that metapgram
# num_samples: Number of samples represented in the metaprogram
# num_patients: Number of patients represented in the metaprogram
# num_primary_samples: Number of samples with primary as disease_timing in the metaprogram
# num_mets_samples: Number of samples with Metastasis as disease_timing in the metaprogram
# num_recurrence_sample: Number of samples with Recurrence as disease_timing in the metaprogram
# proportion_of_samples: Proportion of the total samples represented in the metaprogram
# proportion_of_patients: Proportion of the total patients represented in the metaprogram 
# primary_sample_proportion: Proportion of all primary samples represented
# mets_sample_proportion: Proportion of all metastatic samples represented
# recurrence_sample_proportion: Proportion of all recurrence samples represented
# shannon_entropy: shannon's entropy measuring relative abundance of samples in the metaprogram 
# shannon_entropy_norm: Normalized shannon's entropy 
# eff_number: Effective sample size 
# eff_number_norm: Effective sample size normalized by the total samples in a metaprogram
# greater_pvalue_eff_number_norm: Pvalue of the observed effective sample size > background 
# grater_adj_pvalue_eff_number_norm: Adjused pvalue for effective sample size 
# ORA based metrics ---
# num_sig_genesets: Total number of gene sets identified as significant by ORA
# mean_geneset_specificity: Mean specificity of all gene sets assigned to a metaprogram through ORA
# overall_pvalue_geneset_specificity: pvalue of observed specificity being > or < background 
# overall_adj_pvalue_geneset_specificity: adjusted pvalue for gene set specificity 

# Additional entries in the object will include information about the parameters
# used for metaprogram generation 

# Functions --------------------------------------------------------------------

# helper function for calculating entropy
# copied from https://github.com/AlexsLemonade/ews-nf/blob/main/exploratory-notebooks/04-sample-diversity.Rmd
entropy <- function(sample_n){
  sample_freq <- sample_n / sum(sample_n)
  -sum(sample_freq * log(sample_freq))
}

# calculate stats for all metaprograms 
# returns a dataframe with one row for each metaprogram and one column for each stat
# modified from the function in https://github.com/AlexsLemonade/ews-nf/blob/main/exploratory-notebooks/04-sample-diversity.Rmd
calculate_sample_stats <- function(mp_clusters_df, total_samples, sample_df){
  
  # how many total patients do we have
  total_patients <- length(unique(sample_df$participant_id))
  
  # how many of each disease timing
  disease_timing_totals_df <- table(sample_df$disease_timing) |> 
    as.data.frame() |> 
    tidyr::pivot_wider(
      names_from = c("Var1"),
      values_from = c("Freq")
    )
  
  # df with one row per mp 
  # calculates total spectra, entropy, norm entropy, eff species size
  mp_entropy_df <-  mp_clusters_df <- mp_clusters_df |> 
    dplyr::left_join(sample_df, by = c("sample_id" = "scpca_sample_id")) |> 
    dplyr::group_by(metaprogram, sample_id) |>
    # first get total spectra per sample
    dplyr::summarize(
      num_spectra_per_sample = dplyr::n(), 
      # keep the participant and disease timing information for each sample
      participant_id = dplyr::first(participant_id),
      disease_timing = dplyr::first(disease_timing),
      .groups = "drop_last") |>  
    # now calculate stats for each metaprogram, drop sample grouping 
    dplyr::bind_cols(disease_timing_totals_df) |> 
    dplyr::summarize(
      # totals for each group (spectra, sample, patient, disease timing)
      num_spectra_per_mp = sum(num_spectra_per_sample),
      num_samples = length(sample_id),
      num_patients = length(unique(participant_id)),
      num_primary_samples = sum(disease_timing == "Primary"),
      num_mets_samples = sum(disease_timing == "Metastasis"),
      num_recurrence_samples = sum(disease_timing == "Recurrence"),
      # proportions of total for each category
      proportion_of_samples = num_samples / total_samples,
      proportion_of_patients = num_patients / total_patients,
      primary_sample_proportion = num_primary_samples / unique(Primary),
      mets_sample_proportion = num_mets_samples / unique(Metastasis), 
      recurrence_sample_proportion = num_recurrence_samples / unique(Recurrence), 
      # entropy and neff
      shannon_entropy = entropy(num_spectra_per_sample), # entropy
      shannon_entropy_norm = shannon_entropy / log(num_samples),# normalize by total possible samples in that mp
      eff_number = exp(shannon_entropy), # effective species size or Hill number order of 1
      eff_number_norm = eff_number/num_samples # eff number should be close to the total possible samples in that mp
    )
  
  return(mp_entropy_df)
  
}

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

# helper function for computing metaprogram specificity 
# copied from https://github.com/AlexsLemonade/ews-nf/blob/main/exploratory-notebooks/06-gene-set-uniqueness.Rmd
compute_metaprogram_specificity <- function(ora_results_df, n_metaprograms) {
  
  # since not all metaprograms have gene sets, make sure we track all of them
  mp_names <- sprintf("MP%02d", seq(1, n_metaprograms))
  
  # calculate specificity for each gene set
  # one row per gene set
  geneset_specificity_df <- ora_results_df |>
    dplyr::group_by(ID) |>
    dplyr::summarize(num_mps = dplyr::n_distinct(metaprogram)) |>
    dplyr::mutate(geneset_specificity = 1 - (num_mps - 1) / (n_metaprograms - 1))
  
  # add gene set specificity to metaprograms
  # get one row per metaprogram with mean gene set specificity 
  data.frame(metaprogram = mp_names) |>
    dplyr::left_join(ora_results_df, by = "metaprogram") |>
    dplyr::left_join(geneset_specificity_df, by = "ID") |>
    dplyr::summarize(mean_geneset_specificity = mean(geneset_specificity), .by = "metaprogram")
  
}

# shuffle gene sets assigned to each metaprogram and then re-calculate specificity 
# weight gene sets based on the number of times they appear in the table
# copied from https://github.com/AlexsLemonade/ews-nf/blob/main/exploratory-notebooks/06-gene-set-uniqueness.Rmd
permute_specificity <- function(
    ora_results_df,
    n_metaprograms,
    nreps
){
  
  # get gene set weights 
  geneset_weights <- table(ora_results_df$ID)
  # pool of all significant gene sets across all metaprograms
  all_genesets <- names(geneset_weights)
  
  # get the size (number of sig gene sets) for each metaprogram
  mp_geneset_counts <- ora_results_df |>
    dplyr::group_by(metaprogram) |>
    dplyr::summarise(geneset_count = dplyr::n_distinct(ID), .groups = "drop")
  
  # shuffle genesets assigned to each metaprogram
  background <- replicate(
    nreps, 
    {
      # keep the total number of genesets the same
      shuffled_ora_results <- mp_geneset_counts$geneset_count |> 
        purrr::set_names(mp_geneset_counts$metaprogram) |>
        purrr::map(\(count){
          # sample without replacement from the full pool
          # ensuring no duplicate gene sets within a metaprogram
          sample(all_genesets, count, replace = FALSE, prob = geneset_weights)
        }) |> 
        tibble::enframe(name = "metaprogram", value = "ID") |>
        tidyr::unnest_longer(ID)
      
      # number of MPs each gene set is observed in 
      geneset_prop_df <- shuffled_ora_results |>
        dplyr::group_by(ID) |>
        dplyr::summarize(
          # get the num mps and proportion that each gene set is significant in
          num_mps = dplyr::n_distinct(metaprogram)
        ) |>
        dplyr::mutate(
          geneset_specificity = 1-(num_mps - 1)/(n_metaprograms - 1)
        )
      
      mp_geneset_counts |> 
        dplyr::left_join(shuffled_ora_results, by = c("metaprogram")) |> 
        dplyr::left_join(geneset_prop_df, by = c("ID")) |> 
        dplyr::summarize(
          mean_geneset_specificity = mean(geneset_specificity),
          .by = "metaprogram"
        ) |> 
        dplyr::mutate(
          # get the overall mean across all MPs
          overall_geneset_specificity = mean(mean_geneset_specificity, na.rm = TRUE)
        ) 
      
    }, simplify = FALSE
  ) |>
    dplyr::bind_rows(.id = "replicate")
  
  return(background)
}

# calculate significance for permutation testing
# reports a pvalue for a specified statistic 
calculate_permutation_significance <- function(stat, observed_df, background_df, nreps){
  
  # combine the observed values with the background
  combined_df <- observed_df |> 
    dplyr::select(
      metaprogram,
      obs_value = {{ stat }}
    ) |> 
    dplyr::left_join(background_df, by = c("metaprogram")) |> 
    dplyr::select(-replicate)
  
  # calculate the pvalue for num genesets 
  pvalue_df <- combined_df |> 
    dplyr::mutate(
      # indicate which rows have background > or < obs
      greater_than_obs = {{ stat }} >= obs_value,
      lower_than_obs = {{ stat }} <= obs_value
    ) |> 
    dplyr::group_by(metaprogram) |> 
    dplyr::summarize(
      # calculate p values and retain a column with the observed value
      obs_value = unique(obs_value),
      greater_pvalue = (sum(greater_than_obs) +1)/(nreps + 1),
      lower_pvalue = (sum(lower_than_obs) + 1)/(nreps + 1),
      overall_pvalue = min(greater_pvalue, lower_pvalue) * 2
    ) |> 
    # add adjusted pvalue 
    # depending on the stat will depend on which pvalue we use, either 1 or 2 sided test
    dplyr::mutate(
      greater_adj_pvalue = p.adjust(greater_pvalue, method = "BH"), 
      lower_adj_pvalue = p.adjust(lower_pvalue, method = "BH"),
      overall_adj_pvalue = p.adjust(overall_pvalue, method = "BH")
    )
  
  return(pvalue_df)
  
}

# Read in arguments ------------------------------------------------------------

library(optparse)
library(SingleCellExperiment)

option_list <- list(
  make_option(
    opt_str = c("--metaprograms_file"),
    type = "character",
    default = here::here("ews-nf-out/metaprograms/metaprograms-10_cnmf-5-10-15_filtered/metaprograms.rds"),
    help = "Path to RDS file containing metaprograms"
  ),
  make_option(
    opt_str = c("--sample_metadata_file"),
    type = "character",
    default = here::here("sample-info/sample-metadata.tsv"),
    help = "Path to TSV file with sample metadata"
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
  "Sample metadata file does not exist" = file.exists(opt$sample_metadata_file),
  "Custom gene sets file does not exist" = file.exists(opt$custom_geneset_file),
  "MSigDB gene sets file does not exist" = file.exists(opt$msigdb_geneset_file),
  "Gene set categories file does not exist" = file.exists(opt$category_file),
  "Output file must end in .rds" = stringr::str_detect(opt$output_file, "\\.rds$")
)

# object with metaprogram results
mp_object <- readr::read_rds(opt$metaprograms_file)

# read in sample metadata
sample_df <- readr::read_tsv(opt$sample_metadata_file) |> 
  # keep sample id and then other columns to look at distribution within metaprograms
  dplyr::select(scpca_sample_id, participant_id, disease_timing, tumor_type)

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

# Difference in medians/ Coherence ---------------------------------------------

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

# get the total number of samples/patients to use as the denominator for proportion
# use all samples/patients even if they are not all present in the spectra after filtering
total_samples <- length(mp_object$unique_ids)

# add sample column 
mp_clusters_df <- mp_lookup |>
  dplyr::mutate(sample_id = stringr::str_extract(spectra, "SCPCS[0-9]{6}")) 

# get observed stats
# includes proportion samples, proption patients, proportion for each disease timing 
# shannon entropy and effective sample size
obs_stats_df <- calculate_sample_stats(mp_clusters_df, total_samples, sample_df)

# get the background
background_mp_stats_df <- replicate(opt$nreps, {
  # for each rep, shuffle the spectra
  # keep the distribution of num spectra/ mp the same
  dplyr::mutate(mp_clusters_df, metaprogram = sample(metaprogram)) |>
    calculate_sample_stats(total_samples, sample_df)
}, simplify = FALSE) |>
  dplyr::bind_rows(.id = "replicate")

# pvalue for neff to add in to the final data frame
neff_pvalue_df <- calculate_permutation_significance(eff_number_norm,
                                                     obs_stats_df,
                                                     background_mp_stats_df,
                                                     opt$nreps) |> 
  # select the columns that we want
  # only care about greater than pvalue
  dplyr::select(
    metaprogram,
    "greater_pvalue_eff_number_norm" = "greater_pvalue",
    "greater_adj_pvalue_eff_number_norm" = "greater_adj_pvalue"
  )

# ORA --------------------------------------------------------------------------

# grab hallmarks to use in addition to the custom and Ews specific gene lists 
hallmarks_msigdb_df <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") |> 
  dplyr::select(gs_name, ensembl_gene)
# TODO: add in Gene ontology as a gene set? 

# get genes in the EwS specific gene sets from MSigDB
ef_msigdb_geneset_df <- prep_msigdb_genesets(ef_msigdb_df)

# create a single term2genedf for running ORA with custom, Ews MSig, and hallmarks
all_term_df <- ef_msigdb_geneset_df |>
  dplyr::bind_rows(custom_geneset_df) |> 
  dplyr::bind_rows(hallmarks_msigdb_df) |> 
  tidyr::drop_na()

# extract top genes
top_genes <- extract_top_genes(mp_object$metaprogram_list)

# get ora results
ora_results_df <- run_ora(
  mp_top_list = top_genes,
  gene_universe = mp_object$gene_universe,
  term2gene_df = all_term_df
)

# table of number of gene sets per metaprogram 
num_genesets_df <- get_total_genesets(names(mp_object$metaprogram_list), ora_results_df)

# Gene set specificity ---------------------------------------------------------

# calculate gene set specificity
obs_specificity_df <- compute_metaprogram_specificity(ora_results_df, mp_object$n_metaprograms)

# Calculate background for gene set specificity 
background_specificity_df <- permute_specificity(ora_results_df, mp_object$n_metaprograms, opt$nreps)

# calculate pvalue for specificity 
specificity_pvalue_df <- calculate_permutation_significance(
  mean_geneset_specificity,
  obs_specificity_df,
  background_specificity_df,
  opt$nreps
  ) |> 
  dplyr::select(
    "metaprogram",
    "overall_pvalue_geneset_specificity" = "overall_pvalue",
    "overall_adj_pvalue_geneset_specificity" = "overall_adj_pvalue"
  )

# Combine and export -----------------------------------------------------------

# join all mp metrics 
mp_metrics_df <- max_cor_df |> 
  dplyr::left_join(mp_medians_df, by = "metaprogram") |> 
  # add in all stats
  dplyr::left_join(obs_stats_df, by = "metaprogram") |> 
  # add pvalue for neff
  dplyr::left_join(neff_pvalue_df, by = "metaprogram") |> 
  # add total number of gene sets for ORA 
  dplyr::left_join(num_genesets_df, by = "metaprogram") |> 
  # gene set specificity and pvalue
  dplyr::left_join(obs_specificity_df, by = "metaprogram") |> 
  dplyr::left_join(specificity_pvalue_df, by = "metaprogram")

# create an object for exporting with mp metrics and ora results
metrics_obj <- list(
  mp_metrics_df = mp_metrics_df,
  ora_results_df = ora_results_df,
  background_mp_stats_df = background_mp_stats_df,
  background_specificity_df = background_specificity_df,
  # add some key parameters which will help when we combine these TSV files across different parameters
  n_metaprograms = mp_object$n_metaprograms,
  cnmf_k_range = mp_object$cnmf_k_range,
  filter_spectra = mp_object$filter_spectra,
  orphan_cutoff = mp_object$orphan_cutoff,
  removed_unique_ids = mp_object$removed_unique_ids
)

readr::write_rds(metrics_obj, opt$output_file)

