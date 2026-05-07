# ews-nf

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents** 

- [Description](#description)
- [Running the workflow](#running-the-workflow)
  - [Building the normal reference](#building-the-normal-reference)
  - [Running the main annotation workflow](#running-the-main-annotation-workflow)
  - [Parameters](#parameters)
  - [Profiles](#profiles)
  - [Running the workflow on the Pittsburgh High Throughput Cluster](#running-the-workflow-on-the-pittsburgh-high-throughput-cluster)
  - [Monitoring the workflow with Seqera](#monitoring-the-workflow-with-seqera)
- [Output files](#output-files)
- [Managing software dependencies](#managing-software-dependencies)
  - [Additional software not managed with `renv`](#additional-software-not-managed-with-renv)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

![Overview of ews-nf](./manuscript/figures/pngs/ews-nf_annotation-overview.png)

## Description 

This repository holds a Nextflow workflow used to analyze single-cell RNA-sequencing samples obtained from patients with Ewing sarcoma. 
The main workflow, used to annotate tumor cells, is organized into the following modules: 

* `build-reference`: This module extracts the counts for any cells that should be included in a normal cell reference and saves to a single `SingleCellExperiment` object for use with other modules. 
An additional TSV file containing the library ID, cell barcode, and consensus cell type annotations is also saved.
Any cells with a consensus cell type present in the table of reference cell types (see [`references/reference-celltypes.tsv`](./references/reference-celltypes.tsv)) will be included in the normal reference. 
* `aucell`: [`AUCell`](https://www.bioconductor.org/packages/release/bioc/html/AUCell.html) is used to calculate the AUC for all provided gene sets on each cell in the input `SingleCellExperiment` object. 
Gene sets used are found in `genesets/msigdb-gene-sets.tsv` and `genesets/custom-gene-sets.tsv`. 
Once `AUCell` is run on each object, the AUC results for all cells belonging to the normal reference are used to build a background distribution. 
The mean and standard deviation of the background distribution are used to obtain a z-score for the AUC results from all cells in each sample. 
* `infercnv`: [`inferCNV`](https://github.com/broadinstitute/inferCNV/wiki) is used to infer copy number variations in all non-reference cells. 
All cells that are classified as reference cells using the `build-reference` module are used as a normal reference to run `inferCNV` on each library. 
A TSV file containing the estimated CNV profile for each cell along with a heatmap displaying the output from `inferCNV` are saved for each library. 
In addition, the total CNVs are summed and the mean and standard deviation of the reference cells are used to obtain a z-score for the total CNVs from all cells in each sample. 
This same approach is used to calculate the z-score for all `proportion_scaled_*` columns that represent commonly observed CNVs in Ewing sarcoma (gain of chr1, gain of chr8, gain of ch12, and loss of chr16). 
* `annotate-tumor-cells`: Annotations are assigned to all cells based on the CNV and AUCell z-scores: 
  * Tumor cells are identified as those with a CNV z-score (total or Ewing specific CNV) > 2 (tier2), > 1 (tier1), or > 0 (tier0). 
  * Any cells that are tumor cells are further classified into `high` or `low` based on expression of `EF-activated` (high) or `EF-repressed` (low) gene sets as indicated by `AUCell`. 
  The z-score for each gene set is determined and the mean z-score for all gene sets categorized as either `EF-activated` or `EF-repressed` is calculated. 
  Cells are classified if the mean z-score is greater than the specified `--marker_gene_cutoff`, which has a default of 2. 
* `summary-report`: This module generates a HTML report summarizing the results from the workflow for each sample. 

The tumor cells output from the above workflow are then used as input to a workflow used to annotate metaprograms across all samples. 
That workflow is organized into the following modules: 

* `cnmf`: This module runs [`cNMF`](https://github.com/dylkot/cNMF) on all tumor cells on a range of k values. 
Here, tumor cells are those that are labeled as tumor by the main workflow with a CNV z-score > 0 (tier0).
These cells are labeled as `tumor` in the `annotation_tier0` column of the TSV file output by the main annotation workflow (`annotations/<sample_id>/<library_id>-<sample_id>_annotations.tsv`)
The consensus NMF program for each value of k is included in the output.

## Running the workflow 

In order to run the workflow you will first need to clone the repository. 

Prior to running this workflow, quantified gene expression data must be obtained by running [`scpca-nf`](https://github.com/AlexsLemonade/scpca-nf). 
The processed `SingleCellExperiment` and `AnnData` objects output from that workflow are then used as input to `ews-nf`.

### Building the normal reference

Prior to annotating the cells in individual samples, the `reference` workflow should be used to build the necessary reference files. 
This includes the pooled normal cell reference files and the gene order file to be used with `InferCNV`. 
The pooled reference will include cells from all samples present in the specified `scpca-nf` output directory (`scpca_nf_out` parameter). 

The following command can be used to run the reference part of the workflow: 

```sh
nextflow run main.nf --workflow reference
```

### Running the main annotation workflow

After creating the reference files with the `reference` workflow, the following command can be run from within the cloned repository and will run the workflow locally with all default parameters (see more below):

```sh
nextflow run main.nf
```

### Running the metaprogram identification workflow

After annotating tumor cells using the main workflow, the `nmf_metaprograms` workflow can be used to identify metapgrograms across a group of samples. 
This workflow requires the processed `AnnData` objects (`processed_rna.h5ad` files) output from `scpca-nf` and the tumor cell annotations TSV (`annotations/<sample_id>/<library_id>_annotations.tsv`) file output from the main `ews-nf` workflow. 

The following command can be used to run the NMF part of the workflow: 

```sh
nextflow run main.nf --workflow nmf_metaprograms
```

### Parameters

Below is a table of parameters that can be modified: 

| | | 
| -- | -- | 
| `scpca_nf_out` | This parameter is used to define the path to the directory containing output from `scpca-nf`. The default directory is `scpca-nf-out`. | 
| `sample_ids` | Comma or semi-colon separated list of sample ids (names of folders in `scpca_nf_out` containing the `_processed.rds` files) to run the workflow on. By default, the workflow will be run on all samples that are present in `scpca_nf_out`. |
| `outdir` | Full path to the directory to save the output files from `ews-nf`. The default directory is `ews-nf-out` | 
| `skip_infercnv` | Whether or not to skip running inferCNV for samples with existing inferCNV results in the `outdir`. By default, inferCNV will be run for all samples, unless `--skip_infercnv` is used at the command line | 
| `skip_cnmf` | Whether or not to skip running cNMF for samples with existing rNMF results in the `outdir`. By default, cNMF will be run for all samples, unless `--skip_cnmf` is used at the command line | 

To run the workflow specifying these parameters use the following command: 

```sh
nextflow run main.nf \
  --scpca_nf_out <path to scpca-nf output> \
  --sample_ids <comma separated list of sample ids> \
  --outdir <path to directory to save output from workflow>
```

### Profiles 

`ews-nf` can be run using three different profiles: 

- `standard`: This is the default profile. 
If no profile is specified, then the workflow will use this profile and run the workflow locally.

To use the `standard` profile run: 

```sh
nextflow run main.nf
```

- `pitt_hpc`: This profile can be used to run the workflow using the Pittsburgh High Throughput Cluster (HTC). 
**Note:** It is not recommended to run the workflow from the login node on the HTC. 
Instead, please follow the [instructions on running the workflow using the HTC below](#running-the-workflow-on-the-pittsburgh-high-throughput-cluster). 

To use the `pitt_htc` profile run: 

```sh
nextflow run main.nf -profile pitt_htc
```

- `stub`: This profile can be used for testing the workflow configuration and process input/outputs. 

To use the `stub` profile run: 

```sh
nextflow run main.nf -profile stub -stub 
```

### Running the workflow on the Pittsburgh High Throughput Cluster

The `slurm-submit.sh` script can be used to run `ews-nf` on the Pittsburgh High Throughput Cluster (HTC) using the `pitt_htc` profile. 
This script runs the workflow with `--workflow reference` to build the pooled reference for all samples and then runs the main annotation workflow. 

You will need to clone the repository into your home directory and then use the following command to submit the workflow to the SLURM workload manager: 

```sh
sbatch slurm-submit.sh
```

To modify any of the workflow parameters, specify them at the command line when submitting to SLURM: 

```sh
sbatch slurm-submit.sh --scpca_nf_dir [full path to results folder output from scpca-nf] 
```

If the repository is not located in your home directory, be sure to provide the full path to the repo using the `--project_dir` argument: 

```sh
sbatch slurm-submit.sh --project_dir [full path to cloned repository]
```

The `--sample_ids` argument can be used to specify a subset of samples to annotate: 

```sh
sbatch slurm-submit.sh --sample_ids "SCPCSXXXXXX,SCPCSXXXXXX"
```

To skip running inferCNV for any samples with existing inferCNV results, use: 

```sh
sbatch slurm-submit.sh --skip_infercnv 1
```

By default, if the reference files exist, the reference building (`--workflow reference`) will not be re-run. 
To overwrite an existing reference, use `--overwrite_reference 1`

```sh
sbatch slurm-submit.sh --overwrite_reference 1
```

To run the NMF part of the workflow, use the `slurm-submit-nmf.sh` script after running the main `ews-nf` workflow with `slurm-submit.sh`. 
This script can be run using the following command: 

```sh
sbatch slurm-submit-nmf.sh
```

To skip running cNMF for any samples with existing cNMF results, use: 

```sh
sbatch slurm-submit-nmf.sh --skip_cnmf 1
```


### Monitoring the workflow with Seqera 

To track the progress of the workflow run, you can optionally use the [Seqera platform](https://tower.nf/).
Seqera provides some nice visualizations of the resource usage of each job, as well as a easy to navigate record of run logs that can be useful for debugging.
You will need to create an account on https://tower.nf/ (it is recommended to use your github account for this).

You can then follow the instructions on running [Nextflow using the `-with-tower` flag](https://docs.seqera.io/platform-cloud/getting-started/deployment-options#nextflow--with-tower) to create a new token.
Next, you will want to create an environment variable with this token value using the following command at the command line: `export TOWER_ACCESS_TOKEN=<MYTOWERTOKEN>`
We recommend that you also add this line to your `~/.bash_profile` or equivalent (`~/.zshrc`).

To run the workflow using Seqera, use the following command: 

```sh
nextflow run main.nf -with-tower
```

Any runs of the workflow on the Pittsburgh HTC using the `pitt_htc` profile will appear on the Seqera Platform in the `ews-nf` workspace here: https://cloud.seqera.io/orgs/SeanKarlCohort/workspaces/ews-nf/watch 

You must be a member of the `SeanKarlCohort` organization and have access to the `ews-nf` workspace to monitor runs.

## Output files 

Running the workflow with the `reference` workflow will output three reference files found in `<output_dir>/reference`: 

- A gene order file for use with `inferCNV`. 
The name of the file will be representative of the Ensembl version used to create the file, e.g., `Homo_sapiens.GRCh38.110.gene_order.txt`. 
- A TSV file with the cell types for all cells included in the pooled reference, `all-ref-annotations.tsv`.
- An RDS file with a `SingleCellExperiment` object containing the counts for all cells included in the pooled reference, `all-ref-sce.rds`

The main annotation workflow will include a folder with results from each module organized as follows: 

```tree
ews-nf-out
├── annotations
│   └── <sample_id>
│       └── <library_id>-<sample_id>_annotations.tsv
├── aucell
│   └── <sample_id>
│       └── <library_id>-<sample_id>_aucell-results.tsv
├── infercnv
│   └── <sample_id>
│       ├── <library_id>-<sample_id>_infercnv-heatmap.png
│       └── <library_id>-<sample_id>_infercnv-metadata.tsv
└──  reports
    └── <sample_id>
        └── <library_id>-<sample_id>_summary-report.html

```

The `<library_id>-<sample_id>_summary-report.html` file contains plots and tables to visualize the results from the workflow. 
The summarized results, including the final annotations, are available in the `<library_id>-<sample_id>_annotations.tsv` file, which contains the following columns: 

| column | contents | 
| ------ | -------- | 
| `barcodes` | Cell barcode | 
| `unique_id` | Unique sample identifier in the format, `<library_id>-<sample_id>` | 
| `cell_id` | Unique cell identifier in the format, `<library_id>-<sample_id>-<barcode>` | 
| `consensus_celltype_annotation` | Original consensus cell type annotation assigned by `scpca-nf` |
| `ref_or_query` | Whether the cell is part of the pooled reference (`reference`) or query cell (`query`) | 
| `EF-activated` | Mean z-score of all `EWS::FLI1` activated gene sets output by `AUCell` | 
| `EF-repressed` | Mean z-score of all `EWS::FLI1` repressed gene sets output by `AUCell` | 
| `EMT-ECM` | Mean z-score of all EMT and ECM related gene sets output by `AUCell` | 
| `cnv_zscore` | Total CNV z-score | 
| `chr1_gain_zscore` | Proportion of chr1 with CNV gain z-score | 
| `chr8_gain_zscore` | Proportion of chr8 with CNV gain z-score | 
| `chr12_gain_zscore` | Proportion of chr12 with CNV gain z-score | 
| `chr16_loss_zscore` | Proportion of chr16 with CNV loss z-score | 
| `tumor_class` | Indicates if cells are classified as tumor using CNV > 2 (`tumor_tier2`), CNV > 1 (`tumor_tier1`), or CNV > 0 (`tumor_tier0`) or classified as normal | 
| `expression_class` | Indicates if cells have a `EF-activated` (`high`) or `EF-repressed` (`low`) z-score greater than the marker gene threshold (default: 2). `NA` is used to indicate if the cell does not pass either threshold |
| `annotation_tier0` | Combined annotation based on `tumor_class`, `expression_class`, and `consensus_celltype_annotation. Here tumor cells are those that have CNV > 0 | 
| `annotation_tier1` | Combined annotation based on `tumor_class`, `expression_class`, and `consensus_celltype_annotation. Here tumor cells are those that have CNV > 1 | 
| `annotation_tier2` | Combined annotation based on `tumor_class`, `expression_class`, and `consensus_celltype_annotation. Here tumor cells are those that have CNV > 2 | 


## Managing software dependencies 

All R packages are managed with `renv` and a Dockerfile is used to build an image with R, Bioconductor, and all packages included in `renv`. 
This is the default image used for most processes within the workflow that utilize R.

The following command can be used to build the image: 

```sh
docker buildx build . -t ews-nf --platform linux/amd64
```

Alternatively, the latest version should be available via:

```sh
docker pull ghcr.io/alexslemonade/ews-nf:edge
```

To add any packages to the image outside of `renv`, the `Dockerfile` should be updated to include installation of the desired package.
This image is built and pushed to `ghcr.io/alexslemonade/ews-nf:edge` automatically on every update to `ews-nf` (`main`) in which the `Dockerfile` or `renv.lock` files are changed.

### Additional software not managed with `renv`

The `docker` folder contains any additional images that are built for processes that utilize packages that are not managed within `renv`. 

* The `cnmf` folder contains a `pixi.lock` and `Dockerfile` used to build the `ews-nf/cnmf` Docker image. 
This image is used for running the [`cNMF` package](https://github.com/dylkot/cNMF) within the `ews-nf` workflow. 
This image is built and pushed to `ghcr.io/alexslemonade/ews-nf/cnmf:edge` automatically on every update to `ews-nf` (`main`) in which the `docker/cnmf/Dockerfile` or `docker/cnmf/pixi.lock` files are changed.
