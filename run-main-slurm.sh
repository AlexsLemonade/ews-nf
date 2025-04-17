#!/usr/bin/env bash

#SBATCH --job-name=nextflow-ews-nf
#SBATCH --output=slurm-logs/%j.out
#SBATCH --error=slurm-logs/%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=7-00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc

#TODO: ADD IN SAMPLE ID SPECIFICATIONS # 

#####################################################
# This script is used to run the nextflow ews-nf workflow on a slurm cluster

# Any of the following paramaters can be invoked at the command line:

# workflow_version: Version of the ews-nf workflow to run, default is "main"
# scpca_nf_dir : Path to directory with output from scpca-nf 
# outdir: Path to the output directory where the results will be saved
# tower_pat : Personal access token used to track job submission on Nextflow tower

# To submit this script use the following example:
# sbatch run-main-slurm.sh
# sbatch outdir="ews-nf-out" run-main-slurm.sh
######################################################

set -euo pipefail
module load java/18.0.2

# Define the project directory for nextflow
cd $(dirname "$0")
project_dir=$(pwd)

# set defaults for all arguments
workflow_version="main"
scpca_nf_dir="$project_dir/scpca-nf-out"
outdir="$project_dir/ews-nf-out"
tower_pat=${TOWER_ACCESS_TOKEN:-0}
workdir="$SLURM_SCRATCH/nextflow/work"

# make the workdir if it doesn't exist
mkdir -p $workdir

# set up nextflow tower
export TOWER_ACCESS_TOKEN=$tower_pat

tower_flag="-with-tower"
if [[ $tower_pat == 0 ]]; then
    tower_flag=""
fi


# run ews-nf workflow with specified parameters
nextflow run AlexsLemonade/ews-nf \
  -r $workflow_version \
  -profile "pitt_htc" \
  -work-dir $workdir \
  --scpca_nf_dir "$scpca_nf_dir" \
  --outdir "$outdir" \
  $tower_flag
