#!/usr/bin/env Rscript

# This script is used to make the gene order file necessary to run inferCNV 
# on SCE objects processed through scpca-nf 

# This script was adapted from:
# https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/b1b8dea054d30f081279be9277501bd189a3c4d0/analyses/cell-type-ewings/scripts/cnv-workflow/00-make-gene-order-file.R

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--gtf_file"),
    type = "character",
    help = "Path to gtf file to create gene order file."
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    help = "Path to output gene order file"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# check that file exits
stopifnot("GTF file does not exist" = file.exists(opt$gtf_file))

# read in gtf file
gtf <- rtracklayer::import(opt$gtf_file, feature.type = "gene")

# format gene order file
gtf_df <- gtf |>
  as.data.frame() |>
  dplyr::select(gene_id, seqnames, start, end) |>
  dplyr::mutate(seqnames = glue::glue("chr{seqnames}")) |> 
  # only keep chr1 - 22 and chrX and chrY
  dplyr::filter(grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)$", seqnames))

# export gene order file
readr::write_tsv(gtf_df, opt$output_file, col_names = FALSE)
