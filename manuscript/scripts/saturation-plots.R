#!/usr/bin/env Rscript

library(ggplot2)

# Set default ggplot theme
theme_set(
  theme_bw()
)

# Functions ------------

calculate_shared_mps <- function(df, corr_threshold, jaccard_threshold){
  
  # determine if the full metaprogram is shared with the subset
  # do this for each full metaprogram, so that the max number of shared metaprograms
  # is the number of metaprograms in the full cohort
  df_shared <- df |>
    dplyr::mutate(is_shared = top_correlation >= corr_threshold & jaccard >= jaccard_threshold) |>
    dplyr::group_by(subset_id, num_samples) |>
    dplyr::summarise(
      # avoid double counting situations where the full metaprogram maps to more than one subset metaprogram
      n_shared = dplyr::n_distinct(full_mp[is_shared]),
      n_total  = dplyr::n(),
      .groups  = "drop"
    )
  
  return(df_shared)
}

plot_saturation <- function(df_shared, corr_threshold, jaccard_threshold){
  
  ggplot(df_shared, aes(x = factor(num_samples), y = n_shared)) +
    ggforce::geom_sina(size = 2.5, alpha = 0.85, maxwidth = 0.7) +
    geom_smooth(aes(x = as.numeric(as.factor(num_samples))), method = "loess") +
    labs(
      x = "Number of samples",
      y = "Number of shared metaprograms"
    ) +
    ylim(0,8)
}

# File paths -------------

# subset info file to get sample numbers for each subset
subset_info_file <- here::here("sample-info", "sample-subsets.tsv")


# path to the data directory and associated subset results
data_dir <- here::here("ews-nf-out", "metaprograms")

# comparison metrics output by ews-nf
metrics_files <- list.files(
  data_dir,
  pattern = "*_metaprograms-comparison-metrics.tsv", 
  recursive = TRUE, 
  full.names = TRUE
)

# grab the subset ids from the directory for each file
subset_ids <- dirname(metrics_files) |> 
  stringr::word(-1, sep = "/")

# set the names of the files as the subset ids
names(metrics_files) <- subset_ids

# output pdf files
pdf_dir <- here::here("manuscript", "figures", "pdfs")
full_saturation_pdf <- file.path(pdf_dir, "full-saturation-sina.pdf")
filtered_saturation_pdf <- file.path(pdf_dir, "filtered-saturation-sina.pdf")
saturation_heatmap_pdf <- file.path(pdf_dir, "saturation-heatmap.pdf")

# Read in files ----------

# read in subsets and add a column for how many samples in each subset id 
subset_info_df <- readr::read_tsv(subset_info_file, col_types = list(.default = "c")) |> 
  dplyr::count(subset_id, name = "num_samples")

# read in all metrics files and combine into a single df
metrics_df <- metrics_files |> 
  purrr::map(\(x) 
             readr::read_tsv(x, col_types = list(.default = "d", subset_id = "c", subset_mp = "c", full_mp = "c"))
  ) |> 
  dplyr::bind_rows(.id = "subset_id") |> 
  # add the number of samples per subset
  dplyr::left_join(subset_info_df, by = "subset_id")


# define the thresholds to use
corr_threshold <- 0.75
jaccard_threshold <- 0.6

# named vector mapping full_mp codes to descriptive labels
# TODO: update with real labels once we have confirmed names
mp_labels <- c(
  MP01 = "MP01 EWS::FLI1-classic",
  MP02 = "MP02 noise",
  MP03 = "MP03 EWS::FLI1-low",
  MP04 = "MP04 inflammatory",
  MP05 = "MP05 myogenesis",
  MP06 = "MP06 cycling",
  MP07 = "MP07 unknown"
)

# Saturation plots ---------------

# full saturation plot
calculate_shared_mps(metrics_df, corr_threshold, jaccard_threshold) |> 
  plot_saturation(corr_threshold, jaccard_threshold)

ggsave(full_saturation_pdf, width = 7, height = 7)

# remove the noisy metaprogram
filtered_metrics_df <- metrics_df |>
  dplyr::filter(full_mp != "MP02")

# filtered saturation plot 
calculate_shared_mps(filtered_metrics_df, corr_threshold, jaccard_threshold) |> 
  plot_saturation(corr_threshold, jaccard_threshold)

ggsave(filtered_saturation_pdf, width = 7, height = 7)

# Heatmap ------------------

# calculate the percent of each subset with a shared metaprogram for a given full cohort metaprogram
heatmap_df <- metrics_df |>
  dplyr::mutate(
    # if top correlation wasn't calculated then no genes are shared so this should be 0
    top_correlation = ifelse(is.na(top_correlation), 0, top_correlation),
    is_shared = top_correlation >= corr_threshold,
    jaccard >= jaccard_threshold
  ) |>
  dplyr::summarise(
    subset_shares_mp = any(is_shared),   # TRUE if any subset MP matches this full MP
    .by = c("full_mp", "num_samples", "subset_id")
  ) |>
  dplyr::summarise(
    pct_shared = (sum(subset_shares_mp) / dplyr::n()) * 100,
    .by = c("full_mp", "num_samples")
  ) |> 
  # add some meaningful names to each mp 
  dplyr::mutate(
    full_mp = mp_labels[full_mp]
  )

# Build matrix: rows = full_mp, cols = num_samples
mat <- heatmap_df |>
  tidyr::pivot_wider(names_from = num_samples, values_from = pct_shared) |>
  tibble::column_to_rownames("full_mp") |>
  as.matrix()

# Sort columns numerically
mat <- mat[, order(as.numeric(colnames(mat)))]

# Colour scale: white (0%) → dark blue (100%)
col_fun <- circlize::colorRamp2(
  breaks = c(0, 50, 100),
  colors = c("#F8FAFC", "#6366F1", "#1E1B4B")
)

pdf(saturation_heatmap_pdf, width = 7, height = 5) 

# Draw heatmap
ComplexHeatmap::Heatmap(
  mat,
  name                  = "% subsets\nsharing metaprogram",
  col                   = col_fun,
  cluster_rows          = TRUE,
  cluster_columns       = FALSE, 
  show_row_dend         = FALSE,
  row_names_side        = "left",
  row_names_gp          = grid::gpar(fontsize = 11),
  column_names_side     = "bottom",
  column_names_gp       = grid::gpar(fontsize = 11),
  column_names_rot      = 0,
  column_title          = "Number of samples in subset",
  column_title_side     = "bottom",
  column_title_gp       = grid::gpar(fontsize = 13),
  row_title             = "Full cohort metaprogram",
  row_title_gp          = grid::gpar(fontsize = 13),
  rect_gp               = grid::gpar(col = "white", lwd = 1.5),
  heatmap_legend_param  = list(
    title_gp      = grid::gpar(fontsize = 10),
    labels_gp     = grid::gpar(fontsize = 9),
    legend_height = grid::unit(4, "cm"),
    at            = c(0, 25, 50, 75, 100),
    labels        = c("0%", "25%", "50%", "75%", "100%")
  )
) |> 
  ComplexHeatmap::draw(heatmap_legend_side = "right")

dev.off()
