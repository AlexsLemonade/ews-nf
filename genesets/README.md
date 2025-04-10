# Gene sets

This directory contains information for any gene sets that can be used as input to `AUCell`. 
The `setup` folder contains any scripts used to assist in creating these gene sets. 

## MSigDB gene sets

The `msigdb-gene-sets.tsv` file contains a table of gene sets from [`MSigDB`](https://www.gsea-msigdb.org/gsea/msigdb/) that may be relevant to Ewing sarcoma. 

The following columns can be found in the file: 

| Column name | Description | 
| -- | -- | 
| name | Abbreviated name of the gene set | 
| geneset | Full name of gene set as listed in `MSigDB` | 
| category | Broad category for the specified gene set, e.g., Hallmark (H) or Curated gene sets (C2) | 
| subcategory | Sub category for the specified gene set, e.g., CGP or BP | 

## Custom gene sets 

The `custom-gene-sets.tsv` file contains all genes that belong to a custom gene set that was either sent by an investigator or present in a publication. 

The following custom gene sets are included: 

* `anyaud_ews_targets`: A list of the 78 marker genes defined by [Aynaud _et al._](https://doi.org/10.1016/j.celrep.2020.01.049) to be EWS-FLI1 targets.
Figure 4 shows that expression of these targets is correlated with EWS-FLI1 levels at a single-cell level.
We expect these targets to have increased expression in cells with high EWS-FLI1 activity.

* `wrenn_nt5e_targets`: A list of 28 genes from [Wrenn _et al._](https://doi.org/10.1158/1078-0432.CCR-23-1111) that represent the overlap between the top 217 genes correlated with _NT5E_ expression in patient tumors and the top 200 markers of _NT5E+_ Ewing sarcoma cells _in vitro_.
These genes are shown in Figure 5D and 5E.
We expect these targets to have increased expression in cells with low EWS-FLI1 activity.

* `grohar_ews_repressed`: A list of EWS-FLI1 repressed genes curated by Patrick Grohar. 

* `grohar_ews_activated`: A list of EWS-FLI1 activated genes curated by Patrick Grohar. 

The original gene list sent by Patrick Grohar can be found in `Grohar_BESTLIST_forothers2.xlsx`. 
