# helper to read in a list of files with spectra from cNMF and create a single mtx
build_spectra_mtx <- function(spectra_files){
  
  # read in all spectra and combine into a single matrix
  # rows are genes and each column is a spectra
  # any genes that are not present in a spectra will have NA 
  spectra_mtx <- purrr::imap(
    spectra_files, \(file, id) {
      
      # extract k value to use for labeling columns
      k_value <- stringr::str_extract(file, pattern = "k_\\d+") |> 
        stringr::str_extract("\\d+") |>
        as.integer() 
      
      k_value <- sprintf("k%02d", k_value)
      
      # read in file and create dataframe
      df <- read.table(file) |>
        t() |>
        as.data.frame() |>
        tibble::rownames_to_column("gene_id")
      
      # rename columns to have id, k, and nmf number
      nmf_number <- sprintf("%02d", as.integer(colnames(df)[-1]))
      colnames(df)[-1] <- glue::glue("{id}_{k_value}_CNMF{nmf_number}")
      return(df)
    }
  ) |> 
    purrr::reduce(dplyr::full_join, by = "gene_id") |> 
    tibble::column_to_rownames("gene_id") |> 
    as.matrix()
  
  return(spectra_mtx)
  
}


# helper function to build correlation dataframe with one row per pair of spectra
# duplicates are included to aid in calculation of statistics when 
# choosing which spectra to filter 
build_spectra_cor_df <- function(cor_matrix){
  
  # turn the correlation matrix into a dataframe for plotting
  # one row for each correlation 
  all_cor_df <- as.data.frame(as.table(cor_matrix)) |>
    dplyr::rename(spectra1 = Var1, spectra2 = Var2, correlation = Freq) |>
    # these should be characters and not factors
    dplyr::mutate(spectra1 = as.character(spectra1),
                  spectra2 = as.character(spectra2)) |> 
    # remove correlations between the same spectra
    dplyr::filter(spectra1 != spectra2) |> 
    dplyr::mutate(
      # add columns for sample and k value to help distinguish between/within sample comparisons
      sample1 = stringr::word(spectra1, 1, sep = "_"),
      sample2 = stringr::word(spectra2, 1, sep = "_"),
      pair_type = dplyr::if_else(sample1 == sample2, "within_sample", "between_sample")
    )
  
  return(all_cor_df)
  
}

# helper function for identifying spectra to filter out
# returns a list with a vector of redundant spectra (spectra with high within sample correlation)
# and a vector of orphan spectra (low values of cross sample correlation)
# this is useful for filtering across each individual value of k
filter_spectra <- function(
    all_cor_df,
    redundant_cutoff = 0.9,
    orphan_cutoff = 0.3
) {
  
  # compute cross-sample correlation max for each spectra
  cross_sample_max <- all_cor_df |>
    # only keep correlations between samples
    dplyr::filter(pair_type == "between_sample") |>
    # find the max correlation for each spectra
    dplyr::summarise(
      max_sample_cor = max(correlation, na.rm = TRUE),
      .by = spectra1
    )
  
  # pull out any combinations that result in correlation > 0.9 in a single sample
  high_cor_edges <- all_cor_df |>
    dplyr::filter(correlation > redundant_cutoff, spectra1 != spectra2, pair_type == "within_sample") |> 
    dplyr::select(spectra1, spectra2, sample1, correlation)
  
  # assign groups to spectra that are all correlated to each other
  high_cor_graph <- igraph::graph_from_data_frame(high_cor_edges, directed = FALSE)
  components <- igraph::components(high_cor_graph)
  group_df <- tibble::tibble(
    spectra1 = names(components$membership),
    group = components$membership
  )
  
  # get the top spectra for each group 
  top_spectra <- group_df |>
    dplyr::left_join(cross_sample_max, by = "spectra1") |>
    dplyr::group_by(group) |>
    dplyr::slice_max(order_by = max_sample_cor, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |> 
    dplyr::pull(spectra1)
  
  # identify redundant spectra
  # these are included in the high correlation groups and not part of the top spectra
  redundant_spectra <- high_cor_edges |> 
    dplyr::select(spectra1, spectra2) |> 
    tidyr::pivot_longer(everything()) |> 
    dplyr::filter(!value %in% top_spectra) |> 
    dplyr::pull(value) |> 
    unique() |>
    as.character() #make sure they aren't factors
  
  # pull out any spectra with cross sample max that doesn't meet cutoff
  orphan_spectra <- cross_sample_max |>
    dplyr::filter(max_sample_cor < orphan_cutoff) |>
    dplyr::pull(spectra1) |>
    as.character() # make sure they aren't factors 
  
  return(
    list(
      redundant = redundant_spectra,
      orphan = orphan_spectra
    ))
  
}