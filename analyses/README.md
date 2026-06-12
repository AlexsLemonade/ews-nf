# Analyses

This folder contains scripts and data for analyses outside of the `ews-nf` workflow. 

## Scripts

* `goodspeed-prep.R`: This script is used to prepare the individual processed `SingleCellExperiment` objects from [Goodspeed et al.](https://doi.org/10.1158/1078-0432.ccr-24-2040). 
Before running this script you must download the data as described below. 

* `run-ora.R`: This script can be used to run over-representation analysis on a single set of metaprograms. 

## Data 

* `goodspeed`: This folder contains files obtained from [`GSE277083`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277083)
Data should be organized in the following structure before running the `goodspeed-prep.R` script: 

```
goodspeed
├── cellranger_outs
│   └── GSE277083_raw_feature
│       ├── barcodes.tsv.gz # renamed from  GSE277083_raw_feature_bc_matrix_barcodes.tsv.gz to match expected file names for DropletUtils
│       ├── features.tsv.gz # renamed from  GSE277083_raw_feature_bc_matrix_features.tsv.gz to match expected file names for DropletUtils
│       └── matrix.mtx.gz # renamed from  GSE277083_raw_feature_bc_matrix_matrix.mtx.gz to match expected file names for DropletUtils
├── GSE277083_aggregation.csv.gz
├── GSE277083_cell_metadata_primary.txt.gz
├── GSE277083_counts_matrix_primary.txt.gz
```

  * `cellranger_outs`: This folder contains the raw output from Cell Ranger for all samples in the study.
  All barcodes with the same suffix, -1, -2, -3, etc. are from the same sample. 
  * `GSE277083_aggregation_csv.gz`: This file contains the order in which samples were combined to create the aggregated output from Cell Ranger.  
  * `GSE277083_cell_metadata_primary.txt.gz`: This file contains the cell metadata, including annotations, for all filtered cells from primary samples used in the manuscript. 
  * `processed_sce_files`: This folder contains the output from running `goodspeed-prep.R` which includes a single processed `SingleCellExperiment` object for each of the 7 primary tumor samples. 


## References 

This folder contains any reference files needed to run scripts in the `scripts` folder. 

* `Homo_sapiens.GRCh38.110.mitogenes.txt`: This file was downloaded from `s3://scpca-nf-references/homo_sapiens/ensembl-110/annotation/Homo_sapiens.GRCh38.110.mitogenes.txt`
