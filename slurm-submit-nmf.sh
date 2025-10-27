#!/usr/bin/env bash

#SBATCH --job-name=nextflow-ews-nmf
#SBATCH --output=slurm-logs/%j.out
#SBATCH --error=slurm-logs/%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=7-00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc

#####################################################
# This script is used to run the NMF part of the nextflow ews-nf workflow on a slurm cluster
# It assumes that you have cloned the workflow into your home directory, use the `--project_dir` argument to specify a different location

# Any of the following parameters can be invoked at the command line:

# sample_ids: Comma separated list of sample IDs to run the workflow on. 
#  sample IDs can be just the sample ID (SCPCSXXXXXX) or the unique ID (SCPCLXXXXXX-SCPCSXXXXXX)
#  By default, the workflow will be run on all samples in `scpca_nf_dir`.
# scpca_nf_dir : Path to `results` directory output from scpca-nf 
# outdir: Path to the output directory where the ews-nf results will be saved
#  This should also be where the results from running the main ews-nf workflow live 
#  as the annotations from that workflow are used as input here 
# k_components_lower: Lower bound for number of components to use with cNMF, default: 5
# k_components_upper: Upper bound for number of components to use with cNMF, default: 20
# k_step_size: Step size for creating the range of components to use with cNMF, default: 1
# skip_cnmf: Used to skip running cNMF for samples where cNMF results already exist in the output directory
#  By default, cNMF will be run on all samples regardless of if results exist
#  To skip, use --skip_cnmf 1
# project_dir: Path to the ews-nf repo, default is $HOME/ews-nf

# To submit this script use the following example:
# sbatch slurm-submit.sh
# sbatch slurm-submit.sh --scpca_nf_dir /path/to/scpca-nf-out --outdir /path/to/ews-nf-out
# sbatch slurm-submit.sh --sample_ids "SCPCS000001,SCPCS000002"
# sbatch slurm-submit.sh --skip_cnmf 1
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
skip_cnmf=${skip_cnmf:-0}
k_components_lower=${k_components_lower:-5}
k_components_upper=${k_components_upper:-20}
k_step_size=${k_step_size:-1}

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

# flag for cnmf
if [[ $skip_cnmf -eq 1 ]]; then 
    cnmf_flag="--skip_cnmf"
else 
    cnmf_flag=""
fi

# run ews-nf workflow with specified parameters
nextflow run "$project_dir/main.nf" \
  -profile "pitt_htc" \
  -entry nmf_metaprograms \
  --scpca_nf_dir "$scpca_nf_dir" \
  --outdir "$outdir" \
  --k_components_lower $k_components_lower \
  --k_components_upper $k_components_upper \
  --k_step_size $k_step_size \
  $sample_flag \
  $infercnv_flag \
