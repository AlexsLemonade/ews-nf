#!/usr/bin/env Rscript

# Read in a list of TSV files with individual data frames
# outputs a single TSV file 

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--input_tsv_files"),
    type = "character",
    default = "",
    help = "Comma separated list of TSV files with data frames to combine"
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    default = "", 
    help = "Path to TSV file to save the combined data frame"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# get list of all files
tsv_files <- stringr::str_split_1(opt$input_tsv_files, ",")

# make a single df
df <- tsv_files |> 
  purrr::map(readr::read_tsv) |> 
  dplyr::bind_rows()

# export
readr::write_tsv(df, opt$output_file)

