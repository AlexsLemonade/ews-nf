#!/usr/bin/env Rscript

# This script is used to recreate figures from the individual-metaprogram-evaluation.Rmd
# report to be used in figures for the SKC manuscript

# this includes the correlation matrix of all NMF spectra
# and the heatmap showing results from ORA 

library(ggplot2)

# Set up -----------------------------------------------------------------------

# input files
manuscript_dir <- here::here("manuscript")

data_dir <- file.path(manuscript_dir, "data", "metaprograms")
mp_object_file <- file.path(data_dir, "k-7_metaprograms.rds")
mp_ora_file <- file.path(manuscript_dir, "full-cohort-metaprograms", "k-7_ora-results.tsv")

# color palette
mp_palette_file <- file.path(manuscript_dir, "palettes", "metaprogram-palette.tsv")

# output files 
pdf_dir <- file.path(manuscript_dir, "figures", "pdfs")
cor_heatmap_pdf_file <- file.path(pdf_dir, "metaprogram-cor-heatmap.pdf")

ora_heatmap_pdf_file <- file.path(pdf_dir, "metaprogram-ora-genesets-heatmap.pdf")

set.seed(2026)

# Correlation ------------------------------------------------------------------

# read in metaprogram object
mp_obj <- readr::read_rds(mp_object_file)

# define colors for the metaprograms
mp_colors <- readr::read_tsv(mp_palette_file) |> 
  tibble::deframe()

# metaprogram/cluster annotation
annotation <- ComplexHeatmap::HeatmapAnnotation(
  metaprogram = sprintf("MP%02d", mp_obj$clusters),
  col = list(
    metaprogram = mp_colors
  ))


# pull out correlations and reconstruct dendogram
cor_mtx <- mp_obj$cor_matrix
dist_mtx <- as.dist(1-cor_mtx)
tree <- hclust(dist_mtx, method = "ward.D2")
dend <- as.dendrogram(tree)

# export the heatmap
pdf(cor_heatmap_pdf_file, width = 6, height = 6)

ComplexHeatmap::Heatmap(
  cor_mtx, 
  show_row_names = FALSE, 
  show_column_names = FALSE, 
  cluster_rows = dend,
  cluster_columns = dend,
  top_annotation = annotation,
  name = "correlation"
) |> 
  ComplexHeatmap::draw()

dev.off()

# Gene set Heatmap -------------------------------------------------------------

# read in metrics and grab ora results
ora_results_df <- readr::read_tsv(mp_ora_file)

# convert ora results into a heatmap of pvalues
heatmap_mtx <- ora_results_df |> 
  dplyr::mutate(logp = -log10(p.adjust)) |> 
  # only keep the top 5 gene sets
  dplyr::slice_max(logp, n = 5, by = metaprogram) |>
  dplyr::select(metaprogram, ID, logp) |> 
  tidyr::pivot_wider(
    names_from = metaprogram,
    values_from = logp,
    values_fill = 0  # set missing gene sets to 0                         
  ) |> 
  # remove _ and wrap gene set names
  dplyr::mutate(
    ID = stringr::str_replace_all(ID, "_", " ") |> 
      stringr::str_to_lower() |>
      stringr::str_wrap(width = 45)
  ) |> 
  tibble::column_to_rownames("ID") |> 
  as.matrix()

ora_ht <- ComplexHeatmap::Heatmap(
  heatmap_mtx,
  name = "-log10(pvalue)\n",
  row_names_side = "left",
  row_dend_side = "right",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  col = circlize::colorRamp2(c(0, 40), colors = c("gray95", "darkslateblue")),
  rect_gp = grid::gpar(col = "white", lwd = 2), # add some white lines around each box
  row_names_gp = grid::gpar(fontsize = 16),
  column_names_gp = grid::gpar(fontsize = 16),
  border = TRUE,
  heatmap_legend_param = list(
    title_gp = grid::gpar(fontsize = 16),
    labels_gp = grid::gpar(fontsize = 16),
    at = c(0, 10, 20, 30, 40), # custom ticks and labels to indicate it goes higher than 40
    labels = c("0", "10", "20", "30", "≥40")
  )
)

pdf(ora_heatmap_pdf_file, width = 12, height = 12)

ora_ht |> 
  # don't cut off super long gene set names
  ComplexHeatmap::draw(padding = unit(c(2, 10, 2, 2), "cm"))

dev.off()
