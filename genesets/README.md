# Gene sets

This directory contains information for any gene sets that can be used as input to `AUCell`. 

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
The `setup` folder contains a script, `convert-symbols.R`, that was used to convert the provided gene symbols to Ensembl identifiers. 

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

* `stag2_loss_down` and `stag2_loss_up`: A list of 232 genes differentially expressed in A673 cells with and without STAG2 as characterized in [Gim enez-Llorente _et al._](https://doi.org/10.1038/s44319-024-00303-6). 
Genes included here were listed in Dataset EV1C.
`stag2_loss_down` includes genes that are downregulated with knock out of STAG2 and `stag2_loss_up` includes genes that were upregulated with knock out of STAG2.

## Gene set categories 

The `gene-set-categories.tsv` file assigns an overall category to all gene sets (both from `MSigDB` and custom) that are used with `AUCell` and for over representation analysis of metaprograms. 
The categories are described below: 

* `EF-activated`: These gene sets contain genes that are activated by the fusion gene, `EWS-FLI1`. 
* `EF-repressed`: These gene sets contain genes that are repressed by the fusion gene, `EWS-FLI1`. 
* `EMT-ECM`: These gene sets contain genes that are found to be highly expressed in tumor cells with low expression of `EWS-FLI`.
* `STAG2-loss`: These gene sets are associated with loss of STAG2 in A673 cells.
* `chemo-resistance`: These gene sets are associated with resistance to chemotherapy drugs typically included in standard of care Ewing sarcoma treatment.
 
## EwS multiome gene sets 

The `Furlan_multiome_3_subgroup_links.csv` file contains all genes associated with each module in the [Ewing sarcoma multi-ome paper](https://doi.org/10.1101/2025.06.18.660457).
Here, modules are representative of potential cell states present in Ewing sarcoma cell lines. 
This file includes all genes output for each of the three modules from [`ArchR::plotPeak2GeneHeatmap()`](https://www.archrproject.com/reference/plotPeak2GeneHeatmap.html) with a k value of 3.
Downstream analyses of the modules utilized the unique genes for each module. 
Original code for generating this file can be found here: https://github.com/furlan-lab/EwS_multiome/blob/669d6ef60a7946ece289ecb5ba353e7fa2c8a2c3/code/fig2.Rmd#L115. 

The `multiome-module-genes.tsv` contains all unique genes associated with each module and the associated Ensembl gene id. 
This file was prepared using `setup/prep-multiome-genesets.R`.
