#!/usr/bin/env Rscript

# This script reads in the data files available as part of GSE277083 and prepares them for ews-nf nmf_metaprograms
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277083
# A full description of the files that have been downloaded can be found in the `analyses/README.md`

# The raw cellranger output is read in, split by sample and then each sample is processed: 

# empty droplets are filtered using DropletUtils::emptyDrops() 
# cells with fewer than 200 genes detected are removed (miQC failed)
# normalization is performed
# annotations for tumor cells are added to the objects 

# Each object is then exported as an individual processed SCE 

library(SingleCellExperiment)

set.seed(2026)

# Functions --------------------------------------------------------------------

# function to filter empty droplets from a SCE 
# returns the filtered SCE object with droplets with <= 0.01 FDR retained
filter_droplets_sce <- function(sce){
  
  # filter empty droplets
  # getting errors about calculating ambient profiles with emptyDropsCellRanger so will use emptyDrops instead
  droplet_stats <- DropletUtils::emptyDrops(counts(sce))
  
  cells_to_retain <- which(droplet_stats$FDR <= 0.01)
  
  sce <- sce[, cells_to_retain]
  
  return(sce)
  
}

# adds per cell QC stats and removes low quality cells
# cells with < 200 detected genes are removed
# all samples here failed miQC or only had 1 cell remaining so instead we use the back up plan from scpca-nf 
filter_low_quality_sce <- function(sce, mito_genes){
  
  sce <- scuttle::addPerCellQCMetrics(
    sce,
    subsets = list(mito = mito_genes[mito_genes %in% rownames(sce)]),
  )
  
  # remove cells that have total genes detected < 200
  # miQC fails with all samples or only gives 1 cell back 
  passing_cells <- sce$detected >= 200
  sce <- sce[, passing_cells]
  
  return(sce)
  
}

# function to normalize SCE objects using the deconvolution method 
# one sample has < 100 cells so we include the backup version of just log-normalization 
# this also mirrors the approach in scpca-nf 
normalize_sce <- function(sce){
  qclust <- NULL
  try({
    # try and cluster similar cells
    # clustering will fail if < 100 cells in dataset
    qclust <- scran::quickCluster(sce)
    
    # Compute sum factors for each cell cluster grouping
    sce <- scran::computeSumFactors(sce, clusters = qclust)
    
    # Include note in metadata re: clustering before computing sum factors
    metadata(sce)$normalization <- "deconvolution"
  })
  
  if (is.null(qclust)) {
    # Include note in metadata re: failed clustering
    metadata(sce)$normalization <- "log-normalization"
  }
  
  # Normalize and log transform RNA counts
  sce <- scuttle::logNormCounts(sce)
  
  return(sce)
}

# adds new columns to the colData of an SCE object with the cell type annotations
# include all potentially relevant annotation columns 
add_annotations <- function(sce, id, annotations_df) {
  
  # filter the annotations to just have the sample that we are looking at
  # theoretically barcodes have suffixes so shouldn't overlap, but seems good to be safe here 
  sample_annotations_df <- annotations_df |>
    dplyr::filter(sample_id == id)
  
  # add in the annotations 
  coldata_df <- colData(sce) |> 
    as.data.frame() |> 
    # we want all column names to be lowercase
    # ews-nf expects plural barcodes for cnmf
    dplyr::rename(barcodes = Barcode) |> 
    dplyr::left_join(sample_annotations_df, by = "barcodes") |> 
    dplyr::mutate(
      celltype_annotation = dplyr::if_else(is.na(celltype_annotation), "Not annotated", celltype_annotation),
      broad_celltype_annotation = dplyr::if_else(is.na(broad_celltype_annotation), "Not annotated", broad_celltype_annotation),
      sub_celltype_annotation = dplyr::if_else(is.na(sub_celltype_annotation), "Not annotated", sub_celltype_annotation)
    )
  
  colData(sce) <- DataFrame(coldata_df, row.names = coldata_df$barcode)
  
  return(sce)
  
}

# File paths -------------------------------------------------------------------

# paths to input data files 
data_dir <- here::here("analyses/data/goodspeed")
cellranger_dir <- file.path(data_dir, "cellranger_outs/GSE277083_raw_feature")
cell_metadata_file <- file.path(data_dir, "GSE277083_cell_metadata_primary.txt.gz")
sample_order_file <- file.path(data_dir, "GSE277083_aggregation.csv.gz")

# list of mito genes that correspond to 10X GRCh38 2020-A reference 
# grabbed from s3://scpca-nf-references
mito_file <- here::here("analyses/references/Homo_sapiens.GRCh38.110.mitogenes.txt")

# output directory 
object_results_dir <- file.path(data_dir, "processed_files")

# Prep individual SCE objects --------------------------------------------------

# read in files 
sce <- DropletUtils::read10xCounts(cellranger_dir)
mito_genes <- readLines(mito_file)

# data frame with cell type annotations and other stats associated with final data 
# annotations used in the main figures are in the Cell.Type column 
cell_metadata_df <- readr::read_tsv(cell_metadata_file)

# pull out just the annotation columns we want to add into the objects
annotations_df <- cell_metadata_df |> 
  dplyr::select(
    # rename to match the default barcodes column in the SCE 
    barcodes = "Cell",
    # rename the rest of the columns for legibility 
    sample_id = "Sample",
    celltype_annotation = "Cell.Type",
    broad_celltype_annotation = "Broad.Cell.Type",
    sub_celltype_annotation = "Sub.Cell.Type"
    
  )

# remove any suffixs that are on the sample id
sample_order_df <- readr::read_csv(sample_order_file, col_types = list(.default = "c")) |> 
  dplyr::mutate(
    sample_id = stringr::str_replace(sample_id, "_ReSeqDec2021", "")
  )

# get the list of all sample ids that are kept for the final analysis in the paper 
# these are the only ones that we can identify tumor cells for 
primary_samples <- cell_metadata_df$Sample |>
  unique()
# get the index which corresponds to the suffix for the barcode in the full sce object
sample_idx <- which(sample_order_df$sample_id %in% primary_samples)

# prep output directories and file names 
# construct output directories using sample ids 
sample_dirs <- file.path(object_results_dir, primary_samples)
fs::dir_create(sample_dirs)

# define the full paths to the files
# use sample id for both folder name and file name
# ews-nf can handle this
sce_files <- glue::glue("{sample_dirs}/{primary_samples}_processed.rds")
anndata_files <- glue::glue("{sample_dirs}/{primary_samples}_processed_rna.h5ad")
annotations_files <- glue::glue("{sample_dirs}/{primary_samples}_annotations.tsv")
json_files <- glue::glue("{sample_dirs}/{primary_samples}_metadata.json")

# Process and annotate SCEs ----------------------------------------------------

# extract the sample order from barcode suffix 
coldata_df <- colData(sce) |>
  as.data.frame() |>
  dplyr::mutate(
    sample_number = stringr::word(Barcode, -1, sep = "-")
  )

# add to sce to use for splitting 
sce$sample_number <- coldata_df$sample_number

# split into a list of independent sce objects
sce_list <- sample_idx |>
  as.character() |>
  purrr::map(\(s) sce[, sce$sample_number == s]) |> 
  purrr::set_names(primary_samples)

filtered_sce_list <- sce_list |> 
  purrr::imap(
    \(sce, id){
      
      sce <- filter_droplets_sce(sce) |> # remove empty droplets
        # remove low quality cells (also adds PerCellQC stats)
        filter_low_quality_sce(mito_genes) |> 
        normalize_sce() |>
        # add the cell type annotations
        add_annotations(id, annotations_df)
      
    }
  )

# export processed objects -----------------------------------------------------
purrr::pwalk(
  list(filtered_sce_list, sce_files, anndata_files, annotations_files, json_files),
  \(sce, sce_file, h5ad_file, ann_file, json_file) {

    readr::write_rds(sce, sce_file, compress = "bz2") # export sce
    zellkonverter::writeH5AD(sce, h5ad_file, X_name = "counts", compression = "gzip") # export anndata

    colData(sce) |>
      as.data.frame() |>
      readr::write_tsv(ann_file) # export annotations from coldata

    # export processed cell number to a json file
    # used by ews-nf to check the number of cells prior to processing
    list(
      "processed_cells" = nrow(sce)
    ) |> 
      jsonlite::write_json(json_file, 
                           auto_unbox = TRUE,
                           pretty = TRUE)
    
  }
)
