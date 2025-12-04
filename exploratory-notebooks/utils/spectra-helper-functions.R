# helper functions for working with cNMF spectra

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
# duplicates are included to aid in calculation of statistics
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
