# References

This directory contains information for any references used in the workflow. 

* `reference-celltypes.tsv`: This file contains a table of all consensus cell type annotations and associated Cell Ontology term identifiers to include when building the normal reference. 
Any cells with a consensus cell type present in this table will be included in the normal reference. 

* `consensus-validation-groups.tsv`: This file contains a table of all possible consensus cell type annotations and an associated broad group term (e.g., all T cell subtypes are part of the broader group T cells). 
The `validation_group_annotation` column is used to group cells to create plots shown in the summary report. 
This file was originally copied from the [`cell-type-consensus` module that is part of `OpenScPCA-analysis`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/v0.2.4/analyses/cell-type-consensus/references/consensus-validation-groups.tsv). 

* `reference-validation-groups.txt`: This file contains a list of all validation groups that should be used to build the normal cell type reference file, `reference-celltypes.tsv`. 
Any consensus cell type annotation assigned to the validation groups listed in this file will be a part of the normal reference. 

* `validation-markers.tsv`: This file contains a table of marker genes used for validating consensus cell types. 
This table was copied from the [`cell-type-consensus` module from `OpenScPCA-analysis`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/v0.2.4/analyses/cell-type-consensus/references/validation-markers.tsv) and a full description of how it was created can be found [here](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/cell-type-consensus/references#cell-type-validation-marker-genes). 

The `scripts` directory contains the following scripts used to generate reference files: 

* `create-normal-ref-celltypes-list.R`: This script reads in the `reference-validation-groups.txt` and `consensus-validation-groups.tsv` to create a list of consensus cell type annotations that are used to build the normal reference used with both `AUCell` and `InferCNV`. 
The output of this script is saved to `reference-celltypes.tsv`
