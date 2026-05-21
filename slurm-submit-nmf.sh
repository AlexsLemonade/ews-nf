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
# It assumes that you have cloned the workflow into your home directory, use the `--ews_dir` argument to specify a different location

# Any of the following parameters can be invoked at the command line:

# sample_ids: Comma separated list of sample IDs to run the workflow on. 
#  sample IDs can be just the sample ID (SCPCSXXXXXX) or the unique ID (SCPCLXXXXXX-SCPCSXXXXXX)
#  By default, the workflow will be run on all samples in `scpca_nf_dir`.
# scpca_nf_dir : Path to `results` directory output from scpca-nf 
# outdir: Path to the output directory where the ews-nf results will be saved
#  This should also be where the results from running the main ews-nf workflow live 
#  as the annotations from that workflow are used as input here 
# workdir: Path to the nextflow work directory
# k_components_lower: Lower bound for number of components to use with cNMF, default: 5
# k_components_upper: Upper bound for number of components to use with cNMF, default: 25
# k_step_size: Step size for creating the range of components to use with cNMF, default: 5
# skip_cnmf: Used to skip running cNMF for samples where cNMF results already exist in the output directory
#  By default, cNMF will be run on all samples regardless of if results exist
#  To skip, use --skip_cnmf 1
# cnmf_k_range:  List of k values to use from existing cNMF results when generating metaprograms
#  can be either "all" or a comma separated string, default is "5,10,15"
# filter_spectra: Used to filter orphan spectra prior to generating metaprograms, default is 1 (true)
#  to keep all spectra, use --filter_spectra 0
# n_metaprograms: Comma separated string indicating the range of number of metaprograms to use for finding the optimal k, default is "8,9,10,11"
# Must be at least 4 numbers for proper computation of metaprogram metrics
# subset_mode: Boolean indicating whether to run in subset mode, default is 0 (false)
#  in subset mode, the workflow will run on individual subsets of the data defined by the sample-info/sample-subsets.tsv
# subset_ids: Comma separated string indicating which subsets to run in subset mode, default is empty string which will run all subsets
# score_sces: Boolean indicating whether to score and annotate sces with the optimal value of k, default is 0 (false)
#  to score and annotate sces, use --score_sces 1
# ews_dir: Path to the ews-nf repo, default is $HOME/ews-nf
# skip_metaprogram_generation: Used to skip the nmf_metaprograms workflow, used to run just the shared_metaprograms workflow
#  This will compare the metaprgams across subsets without regenerating the metaprograms


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

# Define the path to ews-nf
ews_dir=$HOME/ews-nf

# path to input and output data 
# default is a directory in the vast node for Cillo lab
home_dir="/vast/acillo/$USER"

# set defaults for all arguments
sample_ids=""
scpca_nf_dir="$home_dir/scpca-nf-out"
outdir="$home_dir/ews-nf-out"
workdir="$home_dir/work"

# cnmf parameters
skip_cnmf=${skip_cnmf:-0}
k_components_lower=${k_components_lower:-5}
k_components_upper=${k_components_upper:-15}
k_step_size=${k_step_size:-5}

# metaprogram parameters
cnmf_k_range=${cnmf_k_range:-"5,10,15"}
filter_spectra=${filter_spectra:-1}
# comma separated list of values
n_metaprograms=${n_metaprograms:-"8,9,10,11"}
score_sces=${score_sces:-0}

# turn on subset mode to run individual subsets instead of the full cohort
subset_mode=${subset_mode:-0}
subset_ids=${subset_ids:-""}
# skip metaprogram generation and just calculate shared metrics 
# only useful if subset_mode is 1
skip_metaprogram_generation=${skip_metaprogram_generation:-0}

# grab arguments from command line
while [ $# -gt 0 ]; do
    if [[ $1 == *'--'* ]]; then
        v="${1/--/}"
        declare $v="$2"
    fi
    shift
done

# validate n_metaprograms has at least 4 values by counting the number of commas
length_n_metaprograms=$(echo "$n_metaprograms" | tr -cd ',' | wc -c )
if [[ $length_n_metaprograms -lt 3 ]]; then
    echo "ERROR: n_metaprograms must contain at least 4 values"
    exit 1
fi

# set flags for running specific samples
if [ -n "$sample_ids" ]; then
    sample_flag="--sample_ids $sample_ids"
else
    sample_flag=""
fi

# set flags for running specific subsets
if [ -n "$subset_ids" ]; then
    subset_id_flag="--subset_ids $subset_ids"
else
    subset_id_flag=""
fi

# flag for cnmf
if [[ $skip_cnmf -eq 1 ]]; then 
    cnmf_flag="--skip_cnmf"
else 
    cnmf_flag=""
fi

# flag for filtering spectra
if [[ $filter_spectra -eq 1 ]]; then 
    filter_spectra_flag="--filter_spectra"
else 
    filter_spectra_flag=""
fi

# flag for subset mode
if [[ $subset_mode -eq 1 ]]; then 
    subset_mode_flag="--subset_mode"
else 
    subset_mode_flag=""
fi

# flag for scoring sces
if [[ $score_sces -eq 1 ]]; then 
    score_sces_flag="--score_sces"
else 
    score_sces_flag=""
fi

if [[ $skip_metaprogram_generation -eq 1 ]]; then 
    # run ews-nf workflow with specified parameters
    nextflow run "$ews_dir/main.nf" \
      -profile "pitt_htc" \
      --workflow nmf_metaprograms \
      -work-dir "$workdir" \
      --scpca_nf_dir "$scpca_nf_dir" \
      --outdir "$outdir" \
      --k_components_lower $k_components_lower \
      --k_components_upper $k_components_upper \
      --k_step_size $k_step_size \
      --cnmf_k_range $cnmf_k_range \
      --n_metaprograms $n_metaprograms \
      $sample_flag \
      $cnmf_flag \
      $filter_spectra_flag \
      $score_sces_flag \
      $subset_mode_flag \
      $subset_id_flag
fi 

if [[ $subset_mode -eq 1 ]]; then 
    # run shared metaprograms workflow
    nextflow run "$ews_dir/main.nf" \
      -profile "pitt_htc" \
      --workflow shared_metaprograms \
      -work-dir "$workdir" \
      --scpca_nf_dir "$scpca_nf_dir" \
      --outdir "$outdir" 
fi
   