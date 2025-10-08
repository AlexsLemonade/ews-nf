#!/usr/bin/env Rscript
# This script runs inferCNV on an SCE object using the provided normal reference
# This script was adapted from: 
# https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/faf86b3330726f78d961445bd9c3818e806a4ba9/analyses/infercnv-consensus-cell-type/scripts/01_run-infercnv.R
#
# Usage:
#
# Rscript 01_run-infercnv.R \
#   --sce_file <input sce file> \
#   --reference_file <input normal reference file> \
#   --output_dir <output directory for infercnv results> \
#   --temp_dir <temp directory to store output files> 
#
# For information on additional arguments, run:
#
# Rscript 02-run-infercnv.R --help
#

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(optparse)
})

option_list <- list(
  make_option(
    opt_str = "--sce_file",
    type = "character",
    help = "Path to the SCE file to run inferCNV on"
  ),
  make_option(
    opt_str = "--reference_file",
    type = "character",
    help = "Path to the normal reference SCE"
  ),
  make_option(
    opt_str = c("--hmm_model"),
    type = "character",
    default = "i6",
    help = "The HMM model to use with inferCNV, either 'i3' or 'i6'. Default is i6."
  ),
  # clustering options. context:
  # https://github.com/broadinstitute/infercnv/wiki/infercnv-tumor-subclusters#tumor-subclustering-by-leiden-clustering-preferred
  make_option(
    opt_str = c("--k_nn"),
    type = "numeric",
    default = 20,
    help = "Number of nearest neighbors to use during leiden clsutering in inferCNV; default is 20"
  ),
  make_option(
    opt_str = c("--leiden_method"),
    type = "character",
    default = "PCA",
    help = "Data on which clustering should be performed in inferCNV, either PCA or simple (expression data itself) for leiden clustering"
  ),
  make_option(
    opt_str = c("--leiden_function"),
    type = "character",
    default = "CPM",
    help = "Objective function to use during leiden clustering"
  ),
  make_option(
    opt_str = c("--leiden_resolution"),
    type = "character", # note this will get converted to numeric later if needed
    default = "auto",
    help = "Resolution value to use with leiden clustering in inferCNV. By default this is automatically chosen, but a numeric value can be provided.."
  ),
  make_option(
    opt_str = c("--gene_order_file"),
    type = "character",
    help = "Path to gene order file as tab delimited .txt file with no column headers.
      Columns are: Ensembl gene id, chr, start, stop."
  ),
  make_option(
    opt_str = "--output_metadata_file",
    type = "character",
    help = "Path to file to save inferCNV results as a TSV"
  ),
  make_option(
    opt_str = "--output_heatmap_file",
    type = "character",
    help = "Path to file to save inferCNV heatmap as a PNG file"
  ),
  make_option(
    opt_str = c("--temp_dir"),
    type = "character",
    help = "Temp directory to store intermediate files when running inferCNV"
  ),
  make_option(
    opt_str = c("-t", "--threads"),
    type = "integer",
    default = 4,
    help = "Number of multiprocessing threads to use."
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = 2025,
    help = "Random seed to set for reproducibility. Note that inferCNV is only reproducible on a given operating system."
  )
)


# Parse options ----------------
opt <- parse_args(OptionParser(option_list = option_list))


# check inputs and outputs 
stopifnot(
  "sce_file does not exist" = file.exists(opt$sce_file),
  "reference_file does not exist" = file.exists(opt$reference_file),
  "gene_order_file does not exist" = file.exists(opt$gene_order_file),
  "output_metadata_file must end in .tsv" = stringr::str_detect(opt$output_metadata_file, "\\.tsv$"),
  "output_heatmap_file must end in .png" = stringr::str_detect(opt$output_heatmap_file, "\\.png$"),
  "hmm_model not properly specified" = opt$hmm_model %in% c("i3", "i6")
)

# make the resolution parameter numeric if it's not auto
resolution <- opt$leiden_resolution
if (resolution != "auto") {
  resolution <- as.numeric(resolution)
  stopifnot(
    "The resolution parameter should either be a number or the string 'auto'." = !is.na(resolution)
  )
}

message("file checks complete")


# Setup files and paths -----------------
set.seed(opt$seed)

# read SCE file
sce <- readr::read_rds(opt$sce_file)

# Read in reference SCE file
ref_sce <- readr::read_rds(opt$reference_file)

# get unique id and construct scratch directory
unique_id <- stringr::word(basename(opt$sce_file), sep = "_", 1)


# define path to output files in scratch directory that will be copied and saved as output from script 
scratch_metadata_file <- file.path(opt$temp_dir, "map_metadata_from_infercnv.txt")
scratch_png <- file.path(opt$temp_dir, "infercnv.png")

message("files read in")

# Prepare input data  ---------------------------------------------

# Update SCE column names to always prefix with unique id for downstream bookkeeping
colnames(sce) <- glue::glue("{unique_id}_{colnames(sce)}")

# Remove any reference cells which also appear in this library
duplicated_cells <- intersect(colnames(ref_sce), colnames(sce))
sce <- sce[, !(colnames(sce) %in% duplicated_cells)]

# Create input matrix by combining library SCE and reference SCE counts assays
raw_counts_matrix <- cbind(
  counts(sce),
  counts(ref_sce)
)

# Check that we removed the duplicates:
stopifnot(
  "Duplicate cells present in input to inferCNV" = !all(duplicated(colnames(raw_counts_matrix)))
)

# Define reference group name for inferCNV
reference_group_name <- "reference"

# Create and export annotation table for inferCNV:
# "unknown" cells are uncharacterized, and "reference" cells are in the reference
annotations_file <- file.path(opt$temp_dir, "infercnv_annotations.txt")
data.frame(cell_id = colnames(raw_counts_matrix)) |>
  dplyr::mutate(annotations = dplyr::if_else(cell_id %in% colnames(sce), "unknown", reference_group_name)) |>
  readr::write_tsv(annotations_file, col_names = FALSE)

message("annotations file created")
# Run infercnv ------------------------

# create the infercnv object
infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = raw_counts_matrix,
  annotations_file = annotations_file,
  delim = "\t",
  gene_order_file = opt$gene_order_file,
  ref_group_name = reference_group_name
)

message("infercnv object created")

# run infercnv
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
  out_dir = opt$temp_dir,
  denoise = TRUE, 
  HMM = TRUE,
  HMM_type = opt$hmm_model, # specifies i3 or i6
  save_rds = FALSE, # don't save the intermediate rds files
  num_threads = opt$threads,
  ###### clustering settings
  leiden_function = opt$leiden_function,
  leiden_method = opt$leiden_method,
  leiden_resolution = resolution
)

message("infercnv run complete, saving results...")

# Save final results -----------------------------------------------------------

# create table with barcodes and CNVs for each chromosome
infercnv::add_to_seurat(
  seurat_obj = NULL,
  infercnv_output_path = opt$temp_dir
)

# add reference information to metadata file and add sum cnv and z-scores
# we have to read in with base R, since there are rownames
cnv_df <- read.table(scratch_metadata_file, header = TRUE, sep = "\t") |>
  # pull out row names into cell_id column
  tibble::rownames_to_column(var = "cell_id") |> 
  # indicator for reference vs query cell based on `subcluster` name
  dplyr::mutate(cell_group = stringr::str_split_i(subcluster, "_", 1))

# calculate the total cnv sum for each cell 
sum_cnv_df <- cnv_df |> 
  tidyr::pivot_longer(
    starts_with("has_cnv_"),
    names_to = "chr",
    values_to = "cnv"
  ) |>
  # sum all values for each
  dplyr::group_by(cell_id) |>
  dplyr::summarize(total_cnv_per_cell = sum(cnv)) |>
  dplyr::ungroup() 

# add back in to larger df with all cnv information and cell types 
all_cnv_df <- cnv_df |> 
  dplyr::left_join(sum_cnv_df, by = "cell_id")

# calculate a z-score for the total CNV and Ewing specific CNVS using the reference cells as mean and sd
# list of all columns to calculate z score for 
# names will eventually be the column names used in the output
cnv_columns <- c(
  "cnv_zscore" = "total_cnv_per_cell",
  "chr1_gain_zscore" = "proportion_scaled_dupli_chr1",
  "chr8_gain_zscore" = "proportion_scaled_dupli_chr8",
  "chr12_gain_zscore" = "proportion_scaled_dupli_chr12",
  "chr16_loss_zscore" = "proportion_scaled_loss_chr16"
)

# make df with cnv zscores 
zscore_df <- cnv_columns |> 
  # for each column calculate the z-score for each cell relative to the reference cells
  purrr::imap(
    \(column, new_column_name){
      
      ref_cnv_sums <- all_cnv_df |> 
        dplyr::filter(cell_group == "reference") |> 
        dplyr::pull({{column}})
      ref_mean <- mean(ref_cnv_sums)
      ref_sd <- sd(ref_cnv_sums)
      
      z_scores <- (all_cnv_df[[column]] - ref_mean)/ref_sd
      
    }
  ) |> 
  # combine into a df with cell id and one column for each zscore (total, chr1, chr8, chr12, and chr16)
  dplyr::bind_cols() |> 
  dplyr::mutate(cell_id = all_cnv_df$cell_id)

# now add in zscore and export
all_cnv_df |> 
  dplyr::left_join(zscore_df, by = "cell_id") |> 
  readr::write_tsv(opt$output_metadata_file)

# copy png file to output directory
fs::file_copy(scratch_png, opt$output_heatmap_file, overwrite = TRUE)
