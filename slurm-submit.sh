#!/usr/bin/env bash

#SBATCH --job-name=nextflow-ews-nf
#SBATCH --output=slurm-logs/%j.out
#SBATCH --error=slurm-logs/%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=7-00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc

#####################################################
# This script is used to run the nextflow ews-nf workflow on a slurm cluster
# It assumes that you have cloned the workflow into your home directory, use the `--project_dir` argument to specify a different location

# Any of the following parameters can be invoked at the command line:

# sample_ids: Comma separated list of sample IDs to run the workflow on. 
#  sample IDs can be just the sample ID (SCPCSXXXXXX) or the unique ID (SCPCLXXXXXX-SCPCSXXXXXX)
#  By default, the workflow will be run on all samples in `scpca_nf_dir`.
# scpca_nf_dir : Path to `results` directory output from scpca-nf 
# outdir: Path to the output directory where the ews-nf results will be saved
# overwrite_reference: Used to re-create the normal reference if files already exist. 
#  The reference will include any samples found in `scpca_nf_dir`. 
# skip_infercnv: Used to skip running inferCNV for samples where inferCNV results already exist in the output directory
#  By default, inferCNV will be run on all samples regardless of if results exist
#  To skip, use --skip_infercnv 1
# project_dir: Path to the ews-nf repo, default is $HOME/ews-nf

# To submit this script use the following example:
# sbatch slurm-submit.sh
# sbatch slurm-submit.sh --scpca_nf_dir /path/to/scpca-nf-out --outdir /path/to/ews-nf-out
# sbatch slurm-submit.sh --sample_ids "SCPCS000001,SCPCS000002"
# sbatch slurm-submit.sh --overwrite_reference 1
######################################################

set -euo pipefail
module load nextflow
# module for using singularity
module load singularity

# Define the project directory for nextflow
project_dir=$HOME/ews-nf

# set defaults for all arguments
sample_ids=""
scpca_nf_dir="$project_dir/scpca-nf-out"
outdir="$project_dir/ews-nf-out"
overwrite_reference=${overwrite_reference:-0}
skip_infercnv=${skip_infercnv:-0}

# grab arguments from command line
while [ $# -gt 0 ]; do
    if [[ $1 == *'--'* ]]; then
        v="${1/--/}"
        declare $v="$2"
    fi
    shift
done

# set flags for running specific samples
if [ -n "$sample_ids" ]; then
    sample_flag="--sample_ids $sample_ids"
else
    sample_flag=""
fi

# flag for infercnv
if [[ $skip_infercnv -eq 1 ]]; then 
    infercnv_flag="--skip_infercnv"
else 
    infercnv_flag=""
fi

# expected reference files 
ref_annotations_file="${outdir}/reference/all-ref-annotations.tsv"
ref_sce_file="${outdir}/reference/all-ref-sce.rds"
gene_order_file="${outdir}/reference/Homo_sapiens.GRCh38.110.gene_order.txt"

# if the ref files don't exist or we are overwriting the reference run the workflow with -entry reference
if [[ ! -f $ref_annotations_file || ! -f $ref_sce_file || ! -f $gene_order_file || $overwrite_reference -eq 1 ]]; then 
    nextflow run "$project_dir/main.nf" \
      -entry reference \
      -profile "pitt_htc" \
      --scpca_nf_dir "$scpca_nf_dir" \
      --outdir "$outdir" \
      --ref_annotations_file "$ref_annotations_file" \
      --ref_sce_file "$ref_sce_file" \
      --gene_order_file "$gene_order_file"

fi

# run ews-nf workflow with specified parameters
nextflow run "$project_dir/main.nf" \
  -profile "pitt_htc" \
  --scpca_nf_dir "$scpca_nf_dir" \
  --outdir "$outdir" \
  $sample_flag \
  $infercnv_flag \
