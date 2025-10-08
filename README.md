# ews-nf

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents** 

- [Description](#description)
- [Running the workflow](#running-the-workflow)
  - [Parameters](#parameters)
  - [Profiles](#profiles)
  - [Running the workflow on the Pittsburgh High Throughput Cluster](#running-the-workflow-on-the-pittsburgh-high-throughput-cluster)
  - [Monitoring the workflow with Seqera](#monitoring-the-workflow-with-seqera)
- [Output files](#output-files)
- [Managing software dependencies](#managing-software-dependencies)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Description 

This repository holds a Nextflow workflow used to analyze single-cell RNA-sequencing samples obtained from patients with Ewing sarcoma. 
The workflow is organized into the following modules: 

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


## Running the workflow 

In order to run the workflow you will first need to clone the repository. 

Prior to running this workflow, quantified gene expression data must be obtained by running [`scpca-nf`](https://github.com/AlexsLemonade/scpca-nf). 
The processed `SingleCellExperiment` objects output from that workflow are then used as input to `ews-nf`. 

### Building the normal reference

Prior to annotating the cells in individual samples, the `reference` entrypoint should be used to build the necessary reference files. 
This includes the pooled normal cell reference files and the gene order file to be used with `InferCNV`. 
The pooled reference will include cells from all samples present in the specified `scpca-nf` output directory (`scpca_nf_out` parameter). 

The following command can be used to run the reference part of the workflow: 

```sh
nextflow run main.nf -entry reference
```

### Running the main annotation workflow

After creating the reference files with the `reference` entrypoint, the following command can be run from within the cloned repository and will run the workflow locally with all default parameters (see more below):

```sh
nextflow run main.nf
```

### Parameters

Below is a table of parameters that can be modified: 

| | | 
| -- | -- | 
| `scpca_nf_out` | This parameter is used to define the path to the directory containing output from `scpca-nf`. The default directory is `scpca-nf-out`. | 
| `sample_ids` | Comma or semi-colon separated list of sample ids (names of folders in `scpca_nf_out` containing the `_processed.rds` files) to run the workflow on. By default, the workflow will be run on all samples that are present in `scpca_nf_out`. |
| `outdir` | Full path to the directory to save the output files from `ews-nf`. The default directory is `ews-nf-out` | 

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
This script runs the workflow with `-entry reference` to build the pooled reference for all samples and then runs the main annotation workflow. 

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

By default, if the reference files exist, the reference building (`-entry reference`) will not be re-run. 
To overwrite an existing reference, use `--overwrite_reference 1`

```sh
sbatch slurm-submit.sh --overwrite_reference 1
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

TBD

## Managing software dependencies 

All software dependencies are managed with `renv` and a Dockerfile is used to build an image with R, Bioconductor, and all packages included in `renv`. 

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
