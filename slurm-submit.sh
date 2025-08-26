#!/usr/bin/env bash

#SBATCH --job-name=nextflow-ews-nf
#SBATCH --output=slurm-logs/%j.out
#SBATCH --error=slurm-logs/%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=7-00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc

#TODO: ADD IN RUNNING SUBSET OF SAMPLE IDS#

#####################################################
# This script is used to run the nextflow ews-nf workflow on a slurm cluster
# It assumes that you have cloned the workflow into your home directory, use the `--project_dir` argument to specify a different location

# Any of the following paramaters can be invoked at the command line:

# scpca_nf_dir : Path to directory with output from scpca-nf 
# outdir: Path to the output directory where the results will be saved
# project_dir: Path to the ews-nf repo, default is $HOME/ews-nf
# tower_pat : Personal access token used to track job submission on Nextflow tower

# To submit this script use the following example:
# sbatch run-main-slurm.sh
# sbatch run-main-slurm.sh --scpca_nf_dir /path/to/scpca-nf-out --outdir /path/to/ews-nf-out
######################################################

set -euo pipefail
module load nextflow
# module for using singularity
module load singularity

# Define the project directory for nextflow
project_dir=$HOME/ews-nf

# set defaults for all arguments
scpca_nf_dir="$project_dir/scpca-nf-out"
outdir="$project_dir/ews-nf-out"
tower_pat=${TOWER_ACCESS_TOKEN:-0}

# grab arguments from command line
while [ $# -gt 0 ]; do
    if [[ $1 == *'--'* ]]; then
        v="${1/--/}"
        declare $v="$2"
    fi
    shift
done

# set up nextflow tower
export TOWER_ACCESS_TOKEN=$tower_pat

tower_flag="-with-tower"
if [[ $tower_pat == 0 ]]; then
    tower_flag=""
fi


# run ews-nf workflow with specified parameters
nextflow run $project_dir/main.nf \
  -profile "pitt_htc" \
  --scpca_nf_dir "$scpca_nf_dir" \
  --outdir "$outdir" \
  $tower_flag
