#!/usr/bin/env Rscript

# This script is used to generate metaprograms using results from cNMF across multiple samples
# Outputs an RDS file with a list of:

# metaprogram_list: List of gene weights for each MP
# clusters: original vector of cluster assignments for each spectra
# cluster_list: Spectra assigned to each cluster/MP
# cor_matrix: Correlation matrix of spectra
# spectra_cor_df: Data frame with all spectra correlation information
# gene_universe: List of all genes used to calculate metaprograms
# n_metaprograms: Number of metaprograms specified
# cnmf_k_range: Set of k values used with cNMF to consider when generating MPs
# filter_spectra: Whether or not orphan spectra were removed
# orphan_cutoff: If orphan spectra were removed, the minimum cross sample correlation
# to keep the spectra. If no filtering, this is -1 
# unique_ids: Vector of unique IDs associated with the MPs

# Functions --------------------------------------------------------------------

# helper to read in a list of files with spectra from cNMF and create a single mtx
build_spectra_mtx <- function(spectra_files){
  
  # read in all spectra and combine into a single matrix
  # rows are genes and each column is a spectra
  # any genes that are not present in a spectra will have NA 
  spectra_mtx <- purrr::imap(
    spectra_files, \(file, id) {
      
      # extract k value to use for labeling columns
      k_value <- stringr::str_extract(file, pattern = "k_\\d+") |> 
        stringr::str_extract("\\d+") |>
        as.integer() 
      
      k_value <- sprintf("k%02d", k_value)
      
      # read in file and create dataframe
      df <- read.table(file) |>
        t() |>
        as.data.frame() |>
        tibble::rownames_to_column("gene_id")
      
      # rename columns to have id, k, and nmf number
      nmf_number <- sprintf("%02d", as.integer(colnames(df)[-1]))
      colnames(df)[-1] <- glue::glue("{id}_{k_value}_CNMF{nmf_number}")
      return(df)
    }
  ) |> 
    purrr::reduce(dplyr::full_join, by = "gene_id") |> 
    tibble::column_to_rownames("gene_id") |> 
    as.matrix()
  
  return(spectra_mtx)
  
}


# helper function to build correlation dataframe with one row per pair of spectra
# duplicates are included to aid in calculation of statistics
build_spectra_cor_df <- function(cor_matrix){
  
  # turn the correlation matrix into a dataframe for plotting
  # one row for each correlation 
  all_cor_df <- as.data.frame(as.table(cor_matrix)) |>
    dplyr::rename(spectra1 = Var1, spectra2 = Var2, correlation = Freq) |>
    # these should be characters and not factors
    dplyr::mutate(spectra1 = as.character(spectra1),
                  spectra2 = as.character(spectra2)) |> 
    # remove correlations between the same spectra
    dplyr::filter(spectra1 != spectra2) |> 
    dplyr::mutate(
      # add columns for sample and k value to help distinguish between/within sample comparisons
      sample1 = stringr::word(spectra1, 1, sep = "_"),
      sample2 = stringr::word(spectra2, 1, sep = "_"),
      pair_type = dplyr::if_else(sample1 == sample2, "within_sample", "between_sample")
    )
  
  return(all_cor_df)
  
}

# generate a character vector of spectra to keep after removing orphans
filter_spectra <- function(cor_matrix, orphan_cutoff){
  
  # create correlation dataframe with info about each spectra
  all_cor_df <- build_spectra_cor_df(cor_matrix)
  
  # compute cross sample max for each spectra
  cross_sample_max <- all_cor_df |>
    # only keep correlations between samples
    dplyr::filter(pair_type == "between_sample") |>
    # find the max correlation for each spectra
    dplyr::summarise(
      max_sample_cor = max(correlation, na.rm = TRUE),
      .by = spectra1
    )
  
  # get a vector of orphan spectra to remove 
  orphan_spectra <- cross_sample_max |>
    dplyr::filter(max_sample_cor < orphan_cutoff) |>
    dplyr::pull(spectra1) |>
    as.character() # make sure they aren't factors 
  
  keep_spectra <- which(!colnames(spectra_mtx) %in% orphan_spectra)
  
  return(keep_spectra)
}

# Read in arguments ------------------------------------------------------------

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--cnmf_results_dir"),
    type = "character",
    # TODO: Remove defaults when adding to nextflow
    default = here::here("ews-nf-out", "cnmf"),
    help = "Path to directory containing spectra files output from cNMF"
  ),
  make_option(
    opt_str = c("--cnmf_k_range"),
    type = "character",
    default = "all",
    help = "Set of k values used with cNMF to keep when clustering spectra into metaprograms, comma separated. Use 'all' to use all available k."
  ),
  make_option(
    opt_str = c("--n_metaprograms"),
    type = "integer",
    default = 10,
    help = "Number of metaprograms to generate"
  ),
  make_option(
    opt_str = c("--filter_spectra"),
    action = "store_true",
    default = FALSE,
    help = "Boolean indicating whether to filter orphan spectra prior to generating metaprograms."
  ),
  make_option(
    opt_str = c("--orphan_cutoff"),
    type = "double",
    default = 0.3,
    help = "Minimum cross sample maximum required to keep spectra if --filter_spectra is used"
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    default = "metaprograms.rds", 
    help = "Path to RDS file to save the object with metaprograms and associated parameters"
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = 2025,
    help = "A random seed for reproducibility."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# set the seed
set.seed(opt$seed)

# make sure inputs are correct
stopifnot(
  "cNMF results directory must be provided" = dir.exists(opt$cnmf_results_dir),
  "orphan cutoff should be between -1 and 1" = dplyr::between(opt$orphan_cutoff, -1, 1),
  "Output file must end in .rds" = stringr::str_detect(opt$output_file, "\\.rds$")
)

if(stringr::str_to_lower(opt$cnmf_k_range) == "all"){
  k_range = c()
} else {
  k_range <- stringr::str_split_1(opt$cnmf_k_range, ",") |> 
    as.integer()
}

# get all the spectra files
spectra_files <- fs::dir_ls(
  path = opt$cnmf_results_dir, 
  glob = "*.gene_spectra_score.k_*.dt_*.txt", 
  recurse = TRUE
) |> 
  # only keep files that have the values of k that we want to keep
  purrr::keep(\(file){
    k <- stringr::str_match(file, "k_(\\d+)")[ ,2] |> 
      as.integer() 
    k %in% k_range || length(k_range) == 0
  })

# extract the unique ID for labeling
unique_ids <- stringr::str_extract(spectra_files, pattern = "SCPCL[0-9]+-SCPCS[0-9]+")
names(spectra_files) <- unique_ids

# check that there are spectra files found and that they have ids 
stopifnot(
  "No spectra files found in the provided cnmf results directory with the specified values of k" = 
    length(spectra_files) != 0,
  "Unable to find unique IDs for each spectra" = length(names(spectra_files)) != 0
)

# Read and filter spectra ------------------------------------------------------

# read in and build spectra mtx
spectra_mtx <- build_spectra_mtx(spectra_files)

# get correlations
cor_mtx <- cor(spectra_mtx, method = "pearson", use = "pairwise.complete.obs")

# filter spectra, identify orphans and remove them
if(opt$filter_spectra){
  
  # list of spectra to keep after removing orphans
  keep_spectra <- filter_spectra(cor_mtx, opt$orphan_cutoff)
  
  # filter spectra mtx and cor mtx
  spectra_mtx <- spectra_mtx[, keep_spectra]
  cor_mtx <- cor_mtx[keep_spectra, keep_spectra]

}

# turn the spectra information into a dataframe 
# save this to the final object for use in the report and calculating metrics
spectra_cor_df <- build_spectra_cor_df(cor_mtx)

# Generate MPs -----------------------------------------------------------------

# cluster spectra into MPs using correlation
dist_mtx <- as.dist(1-cor_mtx)
tree <- hclust(dist_mtx, method = "ward.D2")
dend <- as.dendrogram(tree)
clusters <- cutree(tree, k = opt$n_metaprograms)

# make a list with all spectra in each cluster and name with MP0X
mp_cluster_list <- split(colnames(spectra_mtx), clusters)
names(mp_cluster_list) <- sprintf("MP%02d",seq(1, length(mp_cluster_list)))

# generate list of MPs 
# each MP is a list of gene weights where weights are averaged across all spectra
# names of each list correspond to ensembl gene id
mp_list <- mp_cluster_list |> 
  purrr::map(\(cluster_spectra){
    rowMeans(spectra_mtx[, cluster_spectra, drop = FALSE], na.rm = TRUE)
  })

# Export MPs -------------------------------------------------------------------

# build a list with MPs, clustered spectra, cor mtx and parameters used to export
mp_info_list <- list(
  "metaprogram_list" = mp_list,
  "clusters" = clusters, # original cluster assignments
  "cluster_list" = mp_cluster_list, # spectra grouped by cluster
  "cor_matrix" = cor_mtx,
  "spectra_cor_df" = spectra_cor_df,
  "gene_universe" = rownames(spectra_mtx),
  "n_metaprograms" = opt$n_metaprograms,
  "cnmf_k_range" = opt$cnmf_k_range,
  "filter_spectra" = opt$filter_spectra,
  # if no filtering, this is -1 which means everything was kept
  "orphan_cutoff" = ifelse(opt$filter_spectra, opt$orphan_cutoff, -1),
  "unique_ids" = unique(unique_ids)
)

readr::write_rds(mp_info_list, opt$output_file)
