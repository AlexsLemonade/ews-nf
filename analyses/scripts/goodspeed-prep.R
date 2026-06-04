#!/usr/bin/env Rscript

library(SingleCellExperiment)

set.seed(2026)

# paths to input data files 
data_dir <- here::here("analyses/data/goodspeed")
cellranger_dir <- file.path(data_dir, "cellranger_outs/GSE277083_raw_feature")
cell_metadata_file <- file.path(data_dir, "GSE277083_cell_metadata_primary.txt.gz")
sample_order_file <- file.path(data_dir, "GSE277083_aggregation.csv.gz")

mito_file <- here::here("analyses/references/Homo_sapiens.GRCh38.110.mitogenes.txt")
mito_genes <- readLines(mito_file)

# read in files 
sce <- DropletUtils::read10xCounts(cellranger_dir)
cell_metadata_df <- readr::read_tsv(cell_metadata_file)
# remove any suffixs that are on the sample id
sample_order_df <- readr::read_csv(sample_order_file, col_types = list(.default = "c")) |> 
  dplyr::mutate(
    sample_id = stringr::str_replace(sample_id, "_ReSeqDec2021", "")
  )

# get the list of all sample ids that are kept for the final analysis in the paper 
# these are the only ones that we can identify tumor cells for 
primary_samples <- cell_metadata_df$Sample |>
  unique()
sample_idx <- which(sample_order_df$sample_id %in% primary_samples)

# extract the sample order from barcode suffix 
coldata_df <- colData(sce) |>
  as.data.frame() |>
  dplyr::mutate(
    sample_number = stringr::word(Barcode, -1, sep = "-")
  )

# add to sce 
sce$sample_number <- coldata_df$sample_number

# split into a list of independent sce objects
sce_list <- sample_idx |>
  as.character() |>
  purrr::map(\(s) sce[, sce$sample_number == s]) |> 
  purrr::set_names(primary_samples)


filter_sce <- function(sce){

  # filter empty droplets
  droplet_stats <- DropletUtils::emptyDrops(counts(sce))

  cells_to_retain <- which(droplet_stats$FDR <= 0.01)

  sce <- sce[, cells_to_retain]

  return(sce)

}

sce_list <- sce_list |> 
  purrr::map(filter_sce)


# filter low quality cells using miQC

# filter
