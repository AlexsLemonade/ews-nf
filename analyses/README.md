# Analyses

This folder contains scripts and data for analyses outside of the `ews-nf` workflow. 

## Scripts

## Data 

* `goodspeed`: This folder contains files obtained from [`GSE277083`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277083)

  * `cellranger_outs`: This folder contains the raw output from Cell Ranger for all samples in the study.
  All barcodes with the same suffix, -1, -2, -3, etc. are from the same sample. 
  * `GSE277083_aggregation_csv.gz`: This file contains the order in which samples were combined to create the aggregated output from Cell Ranger.  
  * `GSE277083_cell_metadata_primary.txt.gz`: This file contains the cell metadata, including annotations, for all filtered cells from primary samples used in the manuscript. 
