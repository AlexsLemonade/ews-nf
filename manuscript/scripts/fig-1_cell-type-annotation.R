#!/usr/bin/env Rscript

library(SingleCellExperiment)
library(patchwork)

# Set up -----------------------------------------------------------------------

# input files
manuscript_dir <- here::here("manuscript")

data_dir <- file.path(manuscript_dir, "data")
sample_id <- "SCPCS000845"
sce_file <- file.path(data_dir, sample_id, "SCPCL001228-SCPCS000845_annotated.rds")
celltype_palette_file <-file.path(manuscript_dir, "palettes", "cell-type-palette.tsv")

# output files 
# 1b faceted umaps of AUcell expression and inferCNV totals
png_dir <- file.path(manuscript_dir, "figures", "pngs")
aucell_cnv_png_file <- file.path(png_dir, "1b_infercnv-aucell-umaps.png")

# 1c umap showing cell type annotations
celltype_umap_png_file <- file.path(png_dir, "1c_celltype-umap.png")

# 1d ridge plot showing Ewing specific CNVs in cell types 
cnv_ridgeplot_png_file <- file.path(png_dir, "1d_ewing-cnvs.png")

# Plotting prep ----------------------------------------------------------------

# read in sce
sce <- readr::read_rds(sce_file)

# read in color palette
celltype_colors <- readr::read_tsv(celltype_palette_file, col_types = list(.default = "c")) |> 
  tibble::deframe()

# pull out umap coordinates and specific columns for plotting
umap_df <- scuttle::makePerCellDF(
  sce,
  use.coldata = c(
    "barcodes", 
    "ews_nf_celltype_annotation",
    "cnv_zscore",
    "EF.activated",
    "EF.repressed",
    "chr8_gain_zscore",
    "chr12_gain_zscore"
  ),
  use.dimred = "UMAP"
) |>
  dplyr::rename(
    "UMAP1" = "UMAP.1",
    "UMAP2" = "UMAP.2"
  ) |>
  dplyr::mutate(
    # lump the top cell types since for fig 1c we only show the top 8 cell types 
    ews_lumped_annotation = forcats::fct_lump_n(ews_nf_celltype_annotation, 8, other_level = "other", ties.method = "first") |>
      forcats::fct_infreq() |>
      forcats::fct_relevel("tumor", "tumor-high", "tumor-low", after = 0) |>
      forcats::fct_relevel("unknown", "other", after = Inf),
    # add a tumor category, since we'll use this in 1d
    tumor_category = dplyr::if_else(stringr::str_detect(ews_lumped_annotation, "tumor"), "tumor", "normal")
  )

# Faceted UMAP of AUCell and inferCNV results ----------------------------------

# named vector of columns to plot
# names will be used as facet titles 
columns_to_plot <- c(
  "Total CNV"= "cnv_zscore", 
  "EWS::FLI1 activated" = "EF.activated", 
  "EWS::FLI1 repressed" = "EF.repressed"
)

# faceted umap
panel_b <- columns_to_plot |>
  purrr::imap(\(column, title){

    ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = !!sym(column))) +
      geom_point(alpha = 0.5, size = 0.5) +
      scale_color_viridis_c(option = "turbo") +
      labs(
        color = title
      ) +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1,
        strip.background = element_rect(fill = "transparent"),
        legend.position = "bottom",
        legend.justification = "center",
        legend.title.position = "top"
      )

  }) |>
  patchwork::wrap_plots(ncol = 3) &
  patchwork::plot_annotation(title = sample_id)

ggsave(aucell_cnv_png_file)

# UMAP of final cell type annotations ------------------------------------------

panel_c <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = ews_lumped_annotation)) +
  geom_point(alpha = 0.5, size = 0.2) +
  scale_color_manual(values = celltype_colors) +
  guides(color = guide_legend(
    nrow = 3, 
    override.aes = list(size = 2, alpha = 1)
  )) +
  labs(color = "ews-nf cell type annotation") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    aspect.ratio = 1,
    legend.position = "bottom",
    legend.title.position = "top",
    legend.title = element_text(hjust = 0.5)
  )

ggsave(celltype_umap_png_file)

# Ewing specific CNV plot ------------------------------------------------------

# one row per cell per CNV 
# only show chr8 and chr12 for space 
cnv_df <- umap_df |> 
  tidyr::pivot_longer(cols = c("chr8_gain_zscore", "chr12_gain_zscore"),
                      names_to = "chromosome",
                      values_to = "zscore") |> 
  dplyr::mutate(
    # make the facet labels prettier
    chromosome = stringr::str_remove(chromosome, "_zscore")|> 
      stringr::str_replace("_", " ") |> 
      factor(levels = c("chr8 gain", "chr12 gain"))
  )

# ridge plot showing total CNV for each of the top cell types 
# facet by CNV
panel_d <- ggplot(cnv_df, aes(x = zscore, y = forcats::fct_rev(ews_lumped_annotation), fill = tumor_category)) +
  ggridges::geom_density_ridges(alpha = 0.5, color = "white") +
  scale_fill_manual(values = c("tumor" = celltype_colors[["tumor"]], "normal" = "grey")) +
  facet_wrap(
    vars(chromosome),
    scales = "free_x",
    ncol = 1
  ) +
  labs(
    x = "Proportion CNV",
    y = "",
    fill = ""
  ) +
  theme(legend.position = "top",
        aspect.ratio = 1,
        strip.background = element_rect(fill = "transparent"))

ggsave(cnv_ridgeplot_png_file)
