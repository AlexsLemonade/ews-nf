#!/usr/bin/env Rscript

# helper functions for calculating ORA-based metrics 
# copied from functions in:
# https://github.com/AlexsLemonade/ews-nf/blob/main/modules/metaprograms/resources/usr/bin/03-metaprogram-metrics.R

# function to extract the top genes from each MP 
# returns a list of vectors with the genes for each MP
# grabs top 200 genes by default
extract_top_genes <- function(mp_list, num_genes = 200){
  
  # get the top X number of genes for each MP
  mp_top_list <- mp_list |> 
    purrr::map(\(weights){
      weights |> 
        sort(decreasing = TRUE) |> 
        head(n=num_genes) |> 
        names() # return gene ids not weights
    })
  
  return(mp_top_list)
}

# function for running ORA with clusterProfiler::enricher on each metaprogram
# uses the top 200 genes for each metaprogram by default 
# returns a df with all gene sets for each metaprogram with p.adjust <= 0.05
run_ora <- function(mp_top_list, gene_universe, term2gene_df){
  
  # run ORA using the top 200 genes
  ora_results <- mp_top_list |> 
    purrr::map(\(genes){
      clusterProfiler::enricher(
        gene = genes,
        universe = gene_universe,
        TERM2GENE = term2gene_df
      )
    })
  
  # combine ora results into a single dataframe
  # only include genesets with pvalue < 0.05, the default by enricher 
  ora_results_df <- ora_results |> 
    purrr::map(\(ora) ora@result ) |> 
    purrr::list_rbind(names_to = "metaprogram") |>
    dplyr::filter(p.adjust <= 0.05) 
  
  return(ora_results_df)
  
}

# helper function to calculate the total proportion and number of uniquely mapped gene sets
# geneset_list is a named list of gene sets mapped to each MP after ORA
prop_unique_mapped_genesets <- function(geneset_list){
  
  # count gene sets that appear only once across all metaprograms
  all_observed_genesets <- unlist(geneset_list) 
  total_unique_genesets <- length(unique(all_observed_genesets))
  
  num_unique_genesets <- sum(table(all_observed_genesets) == 1)
  proportion_unique_genesets <- num_unique_genesets/total_unique_genesets
  
  num_geneset_list <- list(
    num_unique_genesets = num_unique_genesets,
    proportion_unique_genesets= proportion_unique_genesets
  )
  
}

# shuffle gene sets assigned to each MP using ORA
# we shuffle geneset assignments ensuring no duplicates per mp
# and calculate total numer and proportion of unique gene sets
# one value per value of k
permute_unique_genesets <- function(ora_results_df, nreps) {
  
  # pool of all significant gene sets across all metaprograms
  all_genesets <- unique(ora_results_df$ID)
  
  # get the size (number of sig gene sets) for each metaprogram
  mp_geneset_counts <- ora_results_df |>
    dplyr::group_by(metaprogram) |>
    dplyr::summarise(geneset_count = dplyr::n_distinct(ID), .groups = "drop")
  
  # shuffle genesets assigned to each metaprogram
  background <- replicate(
    nreps, 
    {
      # keep the total number of genesets the same
      shuffled_groups <- mp_geneset_counts$geneset_count |> 
        purrr::set_names(mp_geneset_counts$metaprogram) |>
        purrr::map(\(count){
          # sample without replacement from the full pool
          # ensuring no duplicate gene sets within a metaprogram
          sample(all_genesets, count, replace = FALSE)
        })
      
      # get unique gene sets as a data frame
      unique_stats <- prop_unique_mapped_genesets(shuffled_groups) |> 
        as.data.frame()
      
    }, simplify = FALSE
  ) |>
    dplyr::bind_rows(.id = "replicate")
  
  return(background)
}
