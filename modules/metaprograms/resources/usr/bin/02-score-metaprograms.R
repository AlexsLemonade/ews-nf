#!/usr/bin/env Rscript

# This script is used to get a score for each metaprogram in each cell in a sample
# the output is a tsv containing the following columns: 
# metaprogram, barcodes, mp_score, unique_id

# Functions --------------------------------------------------------------------

# function for getting a score for a single MP for all cells in an object
score_mp <- function(mp_weights, norm_expr_mtx = norm_expr){
  
  # make sure gene names are present and in the object
  stopifnot(
    "Provided weights are missing gene names" = !is.null(names(mp_weights))
  )
  
  # only consider genes that are present in the sce
  mp_genes <- intersect(names(mp_weights), rownames(norm_expr_mtx))
  mp_weights <- mp_weights[mp_genes]
  
  # check to make sure there are at least some genes present
  stopifnot(
    "No genes found in the SCE object" = length(mp_genes) > 0 
  )
  
  # expression subset: genes × cells
  sub_expr <- norm_expr_mtx[mp_genes, ]
  
  # compute dot product for each cell (matrix multiply)
  scores <- as.numeric(mp_weights %*% sub_expr)
  
  df <- data.frame(
    barcodes = colnames(norm_expr_mtx),
    mp_score = scores
  ) 
  
  return(df)
}

# Read in arguments ------------------------------------------------------------

library(optparse)
library(SingleCellExperiment)

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    default = "",
  ),
  make_option(
    opt_str = c("--annotations_file"),
    type = "character",
    default = "",
    help = "Path to annotations file output from ews-nf with tumor cell annotations"
  ),
  make_option(
    opt_str = c("--unique_id"),
    type = "character",
    default = "",
    help = "Unique Id for library"
  ),
  make_option(
    opt_str = c("--metaprograms_file"),
    type = "character",
    default = "metaprograms.rds",
    help = "Path to RDS file containing metaprograms"
  ),
  make_option(
    opt_str = c("--n_top_genes"),
    type = "integer",
    default = 200, 
    help = "Number of genes to extract for each metaprogram to use for scoring"
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    default = "mp_scores.tsv", 
    help = "Path to TSV file to save cell x score matrix"
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = 2026,
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
  "Processed SCE file does not exist" = file.exists(opt$sce_file),
  "Annotations file does not exist" = file.exists(opt$annotations_file),
  "Metaprograms file does not exist" = file.exists(opt$metaprograms_file),
  "Output file must end in .tsv or .tsv.gz" = stringr::str_detect(opt$output_file, "\\.tsv$|\\.tsv\\.gz$")
)

message("this is a new script")

# read in files
sce <- readr::read_rds(opt$sce_file)
annotations_df <- readr::read_tsv(opt$annotations_file)
mp_obj <- readr::read_rds(opt$metaprograms_file)

# pull out unique Id to save to TSV 
unique_id <- opt$unique_id

# Prep gene expr mtx -----------------------------------------------------------

# get a vector of tumor cells for subsetting
tumor_cells <- annotations_df |> 
  dplyr::filter(stringr::str_detect(annotation_tier0, "tumor")) |> 
  dplyr::pull(barcodes)

# get expression matrix for tumor cells and mp genes of interest 
# only pull out genes that are actually in this object and part of the gene universe
genes_to_use <- intersect(mp_obj$gene_universe, rownames(sce))

# add some temp messages for debugging 
message("SCE dimensions: ", paste(dim(sce), collapse=" x "))
message("genes_to_use length: ", length(genes_to_use))
message("tumor_cells length: ", length(tumor_cells))
message("All genes in SCE: ", all(genes_to_use %in% rownames(sce)))
message("All cells in SCE: ", all(tumor_cells %in% colnames(sce)))
message("logcounts exists: ", "logcounts" %in% assayNames(sce))

expr <- logcounts(sce)[genes_to_use, tumor_cells] 

# quantile normalize the scores 
# preprocessCore requires a regular matrix
norm_expr <- preprocessCore::normalize.quantiles(as.matrix(expr))
# make sure to add back cell and gene info
colnames(norm_expr) <- colnames(expr)
rownames(norm_expr) <- rownames(expr)

# Calculate scores and export --------------------------------------------------

# get the top 200 genes for each mp as a list to use with scoring
mp_top_list <- mp_obj$metaprogram_list |> 
  purrr::map(\(scores){
    scores |> 
      sort(decreasing = TRUE) |> 
      head(n=opt$n_top_genes)
  })

scores_df <- mp_top_list |> 
  purrr::map(score_mp) |> 
  dplyr::bind_rows(.id = "metaprogram") |> 
  dplyr::mutate(unique_id = unique_id)

readr::write_tsv(scores_df, opt$output_file)

