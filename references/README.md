# References

This directory contains information for any references used in the workflow. 

* `reference-celltypes.tsv`: This file contains a table of all consensus cell type annotations and associated Cell Ontology term identifiers to include when building the normal reference. 
Any cells with a consensus cell type present in this table will be included in the normal reference. 
* `consensus-validation-groups.tsv`: This file contains a table of all possible consensus cell type annotations and an associated broad group term (e.g., all T cell subtypes are part of the broader group T cells). 
The `validation_group_annotation` column is used to group cells to create plots shown in the summary report. 
This file was originally copied from the [`cell-type-consensus` module that is part of `OpenScPCA-analysis`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/1e85b284e8533f366de4cba49b50652cf0922539/analyses/cell-type-consensus/references/consensus-validation-groups.tsv). 
* `validation-markers.tsv`: This file contains a table of marker genes used for validating consensus cell types. 
This table was copied from the [`cell-type-consensus` module from `OpenScPCA-analysis`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/1e85b284e8533f366de4cba49b50652cf0922539/analyses/cell-type-consensus/references/validation-markers.tsv) and a full description of how it was created can be found [here](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/cell-type-consensus/references#cell-type-validation-marker-genes). 
