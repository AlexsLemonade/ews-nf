#!/usr/bin/env Rscript

# This script is used to calculate a z-score for the AUC values for each gene set in each cell 
# All AUCell results TSV files are read in along with a table of reference barcodes
# Reference cells are used to define a distribution to use as reference for creating the z-score for all query cells 
# The output is a TSV file for each unique ID with the AUC value and z-score for each cell and each unique gene set
# An additional column is included to indicate if the cell is a reference or query cell 

# Results for all cells provided in the reference will be saved to a separate file 

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--aucell_results_files"),
    type = "character",
    default = "",
    help = "A comma separated list of TSV files containing the results from AUCell"
  ),
  make_option(
    opt_str = c("--reference_cells_file"),
    type = "character",
    help = "Path to TSV with all cells to be included in the reference cells."
  ),
  make_option(
    opt_str = c("--reference_output_file"),
    type = "character",
    help = "Path to TSV file to save the results from AUCell and z-scores for all reference cells"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# make sure input files exist
aucell_files <- stringr::str_split_1(opt$aucell_results_files, ",")
stopifnot(
  "`aucell_results_files` do not exist" = all(file.exists(aucell_files)),
  "Reference cells file must be specified using `--reference_cells_file`" = !is.null(opt$reference_cells_file)
)

# Grab refs and prep AUC df ----------------------------------------------------

# read in reference file and define reference cell ids
ref_cells <- readr::read_tsv(opt$reference_cells_file) |> 
  dplyr::pull(cell_id)

# read in as a single TSV 
aucell_df <- aucell_files |> 
  purrr::map(readr::read_tsv) |>
  dplyr::bind_rows() |> 
  # add ref/query column to aucell 
  dplyr::mutate(cell_id = glue::glue("{unique_id}-{barcodes}"),
                ref_or_query = dplyr::if_else(
                  cell_id %in% ref_cells,
                  "reference",
                  "query"
                ))

# Calculate z-score ------------------------------------------------------------

# get all gene set names 
gene_set_names <- unique(aucell_df$gene_set)

# add zscore based on reference cells for each gene set 
zscore_df <- gene_set_names |> 
  purrr::map(\(name){
    
    # get results for the single gene set
    gene_set_df <- aucell_df |> 
      dplyr::filter(gene_set == name) 
    
    # grab reference cells and get mean and sd
    ref_auc <- gene_set_df |> 
      dplyr::filter(ref_or_query == "reference") |> 
      dplyr::pull(auc)
    
    ref_mean <- mean(ref_auc)
    ref_sd <- sd(ref_auc)
    
    # calculate z-score for each cell 
    gene_set_df <- gene_set_df |> 
      dplyr::mutate(
        auc_zscore = (auc - ref_mean)/ref_sd
      )
    
    return(gene_set_df)
  }) |> 
  dplyr::bind_rows()

# Export results ---------------------------------------------------------------

# save z-scores split by unique_id
zscore_df |> 
  dplyr::group_split(unique_id, .keep = TRUE) |> 
  purrr::walk(\(df) {
    unique_id <- unique(df$unique_id)
    # name output directory by sample ID
    out_dir <- unique(df$sample_id)
    dir.create(out_dir)
    
    out_file <- file.path(out_dir, glue::glue("{unique_id}_aucell-zscores.tsv"))
    readr::write_tsv(df, out_file)
  })

# export all reference cells to a separate file 
zscore_df |> 
  dplyr::filter(ref_or_query == "reference") |> 
  readr::write_tsv(opt$reference_output_file)
