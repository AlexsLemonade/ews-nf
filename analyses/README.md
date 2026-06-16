# Analyses

This folder contains scripts and data for analyses outside of the `ews-nf` workflow. 

## Scripts

* `goodspeed-prep.R`: This script is used to prepare the individual processed `SingleCellExperiment` objects from [Goodspeed et al.](https://doi.org/10.1158/1078-0432.ccr-24-2040). 
Before running this script you must download the data as described below. 

* `run-ora.R`: This script can be used to run over-representation analysis on a single set of metaprograms. 

* `generate-mp-metrics.sh`: This script produces intermediate files needed to run `notebooks/01-SCPCP000015-comparison.Rmd`, which includes calculating metrics between the ScPCA metaprograms and the Sean Karl cohort metaprograms and scoring all ScPCA samples using the Sean Karl cohort metaprogram. 
The script was run with default parameters to generate these files. 

This script also produces intermediate files needed to run `notebooks/02-goodspeed-comparison.Rmd`, which includes calculating metrics between the Goodspeed metaprograms and the Sean Karl Cohort metaprograms and scoring all Goodspeed samples using the Sean Karl cohor metaprograms. 
This command was used to generate these files: 

```sh
project_id="goodspeed" \
  sce_dir="data/goodspeed/processed_files" \
   annotations_dir="data/goodspeed/processed_files" \
   scpca_metaprograms_dir="data/goodspeed/full_cohort" \
   scpca_metaprograms_file="data/goodspeed/full_cohort/k-11_metaprograms.rds" \
   annotations_column="celltype_annotation" \
   celltype_value="Ewing sarcoma" \
   ./SCPCP000015-mp-metrics.sh
```

## Notebooks 

* `01-SCPCP000015-comparison.Rmd`: This notebook compares metaprograms from ScPCA's project SCPCP000015 to the metaprograms from the full Sean Karl cohort. 

* `02-goodspeed-comparison.Rmd`: This notebook compares metaprograms from [`GSE277083`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277083) to the metaprograms from the full Sean Karl cohort.

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
