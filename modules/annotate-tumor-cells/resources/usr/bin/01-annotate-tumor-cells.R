#!/usr/bin/env Rscript

# This script is used to create a table of cell type annotations, including annotation of tumor cells in a single sample

# Tumor cells are first classified by whether or not they have a CNV score or Ewing specific CNV with z-score > 0, 1, or 2
# Tumor cells are further classified by marker gene expression into high and low if they both pass the CNV threshold and have 
# EF-activated > opt$marker_gene_cutoff (high) or EF-repressed > opt$marker_gene_cutoff (low) where the default is 2

# The output contains the following columns: 

# barcodes
# unique_id
# cell_id
# ref_or_query 
# consensus_celltype_annotation

# EF-activated: Mean AUCell z-score for EF-activated gene sets
# EF-repressed: Mean AUCell z-score for EF-repressed gene sets
# EMT-ECM: Mean AUCell z-score for EMT-ECM gene sets
# cnv_zscore: Total CNV zscore 
# chr1_gain_zscore: Chr1 gain zscore 
# chr8_gain_zscore: Chr8 gain zscore
# chr12_gain_zscore: Chr12 gain zscore
# chr16_gain_zscore: Chr16 loss zscore

# tumor_class: Highest CNV threshold passed to be considered tumor, tumor_tier2 means the highest zscore for CNV > 2,
#   tumor_tier1 means the highest zscore for any CNV columns was > 1 but < 2
# expression_class: high (EF-activated) > marker_gene_cutoff or low (EF-repressed) > marker_gene_cutoff
# annotation_tier0: Combined annotation where tumor cells are those that are in all tiers
# annotation_tier1: Combined annotation where tumor cells are those that are in tumor_tier1 or tumor_tier2
# annotation_tier2: Combined annotation where tumor cells are only those that are in tumor_tier2

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--unique_id"),
    type = "character",
    default = "",
    help = "Unique ID"
  ),
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    default = "",
    help = "Path to RDS file containing a processed SingleCellExperiment object."
  ),
  make_option(
    opt_str = c("--aucell_results_file"),
    type = "character",
    default = "",
    help = "Path to TSV file with AUCell results for the given object"
  ),
  make_option(
    opt_str = c("--reference_aucell_file"),
    type = "character",
    default = "",
    help = "Path to TSV file with AUCell reference z-scores"
  ),
  make_option(
    opt_str = c("--gene_set_category_file"),
    type = "character",
    default = "",
    help = "Path to TSV file assigning gene sets used with AUCell to categories object"
  ),
  make_option(
    opt_str = c("--infercnv_results_file"),
    type = "character",
    default = "",
    help = "Path to TSV file with inferCNV metadata table"
  ),
  make_option(
    opt_str = c("--marker_gene_cutoff"),
    type = "numeric",
    default = 2,
    help = "Mean z-score of EF-activated and EF-repressed gene sets to use for classifying cells into tumor-high or tumor-low. 
      Any cells that have a mean z-score greater than this cutoff will be included in the classification."
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    help = "Path to TSV file where annotations will be saved"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# make sure input files exist
stopifnot(
  "sce file does not exist" = file.exists(opt$sce_file),
  "aucell results file does not exist" = file.exists(opt$aucell_results_file), 
  "aucell reference file does not exist" = file.exists(opt$reference_aucell_file),
  "gene set category file does not exist" = file.exists(opt$gene_set_category_file),
  "infercnv results file does not exist" = file.exists(opt$infercnv_results_file)
)

# load SCE
suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

# read in sce 
sce <- readr::read_rds(opt$sce_file)

# read in geneset info 
geneset_df <- readr::read_tsv(opt$gene_set_category_file)

# read in results
aucell_df <- readr::read_tsv(opt$aucell_results_file)
aucell_ref_df <- readr::read_tsv(opt$reference_aucell_file)
infercnv_df <- readr::read_tsv(opt$infercnv_results_file) |> 
  dplyr::select(cell_id, contains("zscore"))

# Get consensus cell types -----------------------------------------------------

# grab coldata with current annotations
coldata_df <- scuttle::makePerCellDF(
  sce,
  use.coldata = c("barcodes", "consensus_celltype_annotation"),
  use.dimred = FALSE
) |> 
  dplyr::mutate(
    unique_id = opt$unique_id,
    cell_id = glue::glue("{unique_id}_{barcodes}")
  ) 

# Combine AUCell results -------------------------------------------------------

# combine all aucell results and get mean for EF-activated, EF-repressed, and EMT-ECM genesets
all_aucell_df <- aucell_df |> 
  # add zscores for reference cells, we want to include all reference cells used to keep AUC and CNV information in one file 
  dplyr::full_join(aucell_ref_df) |> 
  # add gene set categories
  dplyr::left_join(geneset_df, by = c("gene_set" = "geneset_name")) |> 
  # group by geneset category and keep ref_or_query column for output later 
  dplyr::group_by(cell_id, geneset_category, ref_or_query) |> 
  dplyr::summarize(
    mean_category_zscore = mean(auc_zscore)
  ) |> 
  tidyr::pivot_wider(names_from = geneset_category, values_from = mean_category_zscore)

# Annotate ---------------------------------------------------------------------

# combine celltypes, aucell, and infercnv information together
# only keep cells that are in the object, we don't need to classify the reference cells from other libraries
results_df <- all_aucell_df |> 
  dplyr::filter(stringr::str_detect(cell_id, opt$unique_id)) |> 
  # aucell results
  # again using full join here since we want to keep both library unknown and reference cells
  dplyr::left_join(coldata_df, by = "cell_id") |> 
  # infercnv results
  dplyr::left_join(infercnv_df, by = "cell_id") |> 
  dplyr::mutate(
    # first classify the tumor cells by the highest threshold that gets passed across total CNV or any Ewing specific chromosomes
    tumor_class = dplyr::case_when(
      ref_or_query != "reference" & dplyr::if_any(dplyr::ends_with("_zscore"), \(x) x > 2) ~ "tumor_tier2",
      ref_or_query != "reference" & dplyr::if_any(dplyr::ends_with("_zscore"), ~ .x > 1) ~ "tumor_tier1",
      ref_or_query != "reference" & dplyr::if_any(dplyr::ends_with("_zscore"), ~ .x > 0) ~ "tumor_tier0",
      .default = "normal"
    ),
    # column indicating if the cell expresses activated (high) or repressed (low) targets
    expression_class = dplyr::case_when(
      `EF-activated` > opt$marker_gene_cutoff ~ "high",
      `EF-repressed` > opt$marker_gene_cutoff ~ "low",
      .default = NA
    ),
    # combine the annotation to get one column per final annotation per tier
    # if the cell is a tumor cell, note the expression class, otherwise include the consensus cell type annotation 
    annotation_tier0 = dplyr::case_when(
      # any tumor cells should be included here, regardless of class 
      stringr::str_detect(tumor_class, "tumor") & !is.na(expression_class) ~ glue::glue("tumor-{expression_class}"),
      stringr::str_detect(tumor_class, "tumor") ~ "tumor",
      .default = consensus_celltype_annotation
    ),
    annotation_tier1 = dplyr::case_when(
      # here tumor cells are from tier1 and tier2
      (tumor_class == "tumor_tier1" | tumor_class == "tumor_tier2") & !is.na(expression_class) ~ glue::glue("tumor-{expression_class}"),
      tumor_class == "tumor_tier1" | tumor_class == "tumor_tier2" ~ "tumor",
      .default = consensus_celltype_annotation
    ),
    annotation_tier2 = dplyr::case_when(
      # only tier2
      tumor_class == "tumor_tier2" & !is.na(expression_class) ~ glue::glue("tumor-{expression_class}"),
      tumor_class == "tumor_tier2" ~ "tumor",
      .default = consensus_celltype_annotation
    )
  )

# export results 
readr::write_tsv(results_df, opt$output_file)
