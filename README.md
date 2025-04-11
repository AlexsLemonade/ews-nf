# ews-nf

## Description 

This repository holds a Nextflow workflow used to analyze single-cell RNA-sequencing samples obtained from patients with Ewing sarcoma. 
More details on included steps coming soon. 

## Running the workflow 

Prior to running this worklow, quantified gene expression data must be obtained by running [`scpca-nf`](https://github.com/AlexsLemonade/scpca-nf). 
The processed `SingleCellExperiment` objects output from that workflow are then used as input to `ews-nf`. 

The following command will run the workflow locally with all default parameters (see more below):

```sh
nextflow run AlexsLemonade/ews-nf
```

### Parameters

Below is a table of parameters that can be modified: 

| | | 
| -- | -- | 
| `scpca_nf_out` | This parameter is used to define the path to the directory containing output from `scpca-nf`. The default directory is `scpca-nf-out`. | 
| `--sample_ids` | Comma or semi-colon separated list of sample ids (names of folders in `scpca_nf_out` containing the `_processed.rds` files) to run the workflow on. By default, the workflow will be run on all samples that are present in `scpca_nf_out`. |
| `--outdir` | Full path to the directory to save the output files from `ews-nf`. The default directory is `ews-nf-out` | 

To run the workflow specifying these parameters use the following command: 

```sh
nextflow run AlexsLemonade/ews-nf \
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
nextflow run AlexsLemonade/ews-nf 
```

- `pitt_hpc`: This profile can be used to run the workflow using the Pittsburgh HPC. 
**Note that this workflow is still in the process of being set up.**

To use the `pitt_hpc` profile run: 

```sh
nextflow run AlexsLemonade/ews-nf -profile pitt_hpc
```

- `stub`: This profile can be used for testing the workflow configuration and process input/outputs. 

To use the `stub` profile run: 

```sh
nextflow run AlexsLemonade/ews-nf -profile stub -stub 
```

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
