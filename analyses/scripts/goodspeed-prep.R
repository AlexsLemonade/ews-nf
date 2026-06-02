#!/usr/bin/env Rscript

library(SingleCellExperiment)

set.seed(2026)

# paths to input data files 
data_dir <- here::here("analyses/data/goodspeed")
cellranger_dir <- file.path(data_dir, "cellranger_outs/GSE277083_raw_feature")
cell_metadata_file <- file.path(data_dir, "GSE277083_cell_metadata_primary.txt.gz")
sample_order_file <- file.path(data_dir, "GSE277083_aggregation.csv.gz")

# read in files 
sce <- DropletUtils::read10xCounts(data_dir)
sample_order_df <- readr::read_csv(sample_order_file)
cell_metadata_df <- readr::read_tsv(cell_metadata_file)

# extract the sample order from barcode suffix 
coldata_df <- colData(sce) |>
  as.data.frame() |>
  dplyr::mutate(
    sample_number = stringr::word(Barcode, -1, sep = "-")
  )
# add to sce 
sce$sample_number <- coldata_df$sample_number

primary_samples <- c(1, 2, 3, 4, 5, 8, 9, 10, 11, 19)


# split into a list of indepenent samples
sce_list <- primary_samples |>
  as.character() |>
  purrr::set_names() |>
  purrr::map(\(s) sce[, sce$sample_number == s])


filter_sce <- function(sce){

  # filter empty droplets
  droplet_stats <- DropletUtils::emptyDropsCellRanger(
    counts(sce),
    BPPARAM = BiocParallel::MulticoreParam(4))

  cells_to_retain <- which(droplet_stats$FDR <= 0.01)

  sce <- sce[, cells_to_retain]

  return(sce)

}

# filter low quality cells using miQC

# filter
