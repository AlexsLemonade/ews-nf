#!/usr/bin/env Rscript

# This script is used to generate metaprograms using results from cNMF across multiple samples
# Outputs an RDS file with a list of:

# metaprogram_list: List of gene weights for each MP
# cluster_list: Spectra assigned to each cluster/MP
# cor_matrix: Correlation matrix of spectra
# n_metaprograms: Number of metaprograms specified
# cnmf_k_range: Range of k values used with cNMF to consider when generating MPs
# filter_spectra: Whether or not orphan spectra were removed
# orphan_cutoff: If orphan spectra were removed, the minimum cross sample correlation
# to keep the spectra. If no filtering, this is -1 
# unique_ids: Vector of unique IDs associated with the MPs

# source in some helper functions we'll use to prep data 
# build_spectra_mtx() build_spectra_cor_df()
source(here::here("exploratory-notebooks/utils/spectra-helper-functions.R"))

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--cnmf_results_dir"),
    type = "character",
    default = here::here("ews-nf-out", "cnmf"),
    help = "Path to directory containing spectra files output from cNMF"
  ),
  make_option(
    opt_str = c("--cnmf_k_range"),
    type = "integer",
    default = seq(5, 25, 5),
    help = "Path to RDS file containing GeneSetCollection object to use with AUCell."
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
    k %in% opt$cnmf_k_range
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
  
  # create correlation dataframe with info about each spectra
  all_cor_df <- build_spectra_cor_df(cor_mtx)
  
  # compute cross sample max for each spectra
  cross_sample_max <- all_cor_df |>
    # only keep correlations between samples
    dplyr::filter(pair_type == "between_sample") |>
    # find the max correlation for each spectra
    dplyr::summarise(
      max_sample_cor = max(correlation, na.rm = TRUE),
      .by = spectra1
    )
  
  # get a vector of orhpan spectra to keep 
  orphan_spectra <- cross_sample_max |>
    dplyr::filter(max_sample_cor < opt$orphan_cutoff) |>
    dplyr::pull(spectra1) |>
    as.character() # make sure they aren't factors 
  
  keep_spectra <- which(!colnames(spectra_mtx) %in% orphan_spectra)
  
  # filter spectra mtx and cor mtx
  spectra_mtx <- spectra_mtx[, keep_spectra]
  cor_mtx <- cor_mtx[keep_spectra, keep_spectra]

}

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
  purrr::map(\(spectra_list){
    
    rowMeans(spectra_mtx[, spectra_list, drop = FALSE], na.rm = TRUE)
    
  })

# Export MPs -------------------------------------------------------------------

# build a list with MPs, clustered spectra, cor mtx and parameters used to export
mp_info_list <- list(
  "metaprogram_list" = mp_list,
  "cluster_list" = mp_cluster_list,
  "cor_matrix" = cor_mtx,
  "n_metaprograms" = opt$n_metaprograms,
  "cnmf_k_range" = opt$cnmf_k_range,
  "filter_spectra" = opt$filter_spectra,
  # if no filtering, this is -1 which means everything was kept
  "orphan_cutoff" = ifelse(opt$filter_spectra, opt$orphan_cutoff, -1),
  "unique_ids" = unique_ids
)

readr::write_rds(mp_info_list, opt$output_file)
