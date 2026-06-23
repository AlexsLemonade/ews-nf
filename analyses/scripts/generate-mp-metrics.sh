#!/bin/bash

# this script is used to generate files needed for input to analyses/notebooks/01-SCPCP000015-comparison

# This includes:
# 1) A comparison of the metaprograms identified in SCPCP000015 to those identified in the full cohort (using the same metaprogram generation method)
# 2) Scores for each library in SCPCP000015 using the metaprograms identified in the full cohort

# Input arguments: 
# project_id: the ID of the project to analyze
# sce_dir: directory containing processed SCE objects for each library in the specified project
# annotations_dir: directory containing annotation files for each library in the specified project
# cohort_metaprograms_dir: directory containing metaprograms identified in the cohort being compared to SKC
# cohort_metaprograms_file: file containing metaprograms identified in the cohort being compared to SKC
# skc_metaprograms_file: file containing metaprograms identified in the SKC
# annotations_column: column in the annotations file to use for subsetting the data when scoring metaprograms
# celltype_value: value in the annotations_column to use for subsetting the data when scoring metaprograms (e.g., tumor)

# Example usage: 
# project_id="goodspeed" \
#   sce_dir="data/goodspeed/processed_files" \
#   annotations_dir="data/goodspeed/processed_files" \
#   cohort_metaprograms_dir="data/goodspeed/full_cohort" \
#   cohort_metaprograms_file="data/goodspeed/full_cohort/k-11_metaprograms.rds" \
#   annotations_column="celltype_annotation" \
#   celltype_value="Ewing sarcoma" \
#   ./generate-mp-metrics.sh

set -euo pipefail

# navigate to analyses directory to ensure we use the renv file
script_dir=$(dirname "$0")
cd $script_dir/..

# project id
project_id=${project_id:-"SCPCP000015"}

# define paths 
sce_dir=${sce_dir:-"../scpca-nf-out/${project_id}"}
annotations_dir=${annotations_dir:-"data/${project_id}/annotations"}
cohort_metaprograms_dir=${cohort_metaprograms_dir:-"data/${project_id}/metaprograms"}
cohort_metaprograms_file=${cohort_metaprograms_file:-"${cohort_metaprograms_dir}/full_cohort/metaprograms-6_cnmf-5-10-15_filtered/k-6_metaprograms.rds"}
skc_metaprograms_file=${skc_metaprograms_file:-"../ews-nf-out/metaprograms/metaprograms-7_cnmf-5-10-15_filtered/metaprograms.rds"}

# define parameters for scoring metaprograms
annotations_column=${annotations_column:-"annotation_tier0"}
celltype_value=${celltype_value:-"tumor"}

# first compare metrics between the specified metaprograms and SKC metaprograms
Rscript ../modules/compare-metaprograms/resources/usr/bin/01-compare-metaprograms.R \
 --subset_id $project_id \
 --subset_metaprograms_file $cohort_metaprograms_file \
 --full_metaprograms_file $skc_metaprograms_file \
 --output_file ${cohort_metaprograms_dir}/SKC-metaprogram-comparison-metrics.tsv

# now calculate scores for each library using the metaprograms from the full cohort 
sample_ids=$(basename -a ${sce_dir}/*)
for sample_id in $sample_ids; do

    for sce_file in $sce_dir/$sample_id/*_processed.rds; do

        # define library ID
        library_id=$(basename $sce_file | sed 's/_processed.rds$//')

        # grab corresponding annotation file
        annotation_file="${annotations_dir}/${sample_id}/${library_id}-${sample_id}_annotations.tsv"

        # only score if the annotations file exists
        if [ -f $annotation_file ]; then

            # define output file 
            scores_output_file="${cohort_metaprograms_dir}/scores/${library_id}_metaprogram-scores.tsv"

            # generate scores for each library using the full cohort metaprograms
            echo "Generating scores for $library_id"
            Rscript ../modules/metaprograms/resources/usr/bin/02-score-metaprograms.R \
              --sce_file $sce_file \
              --annotations_file $annotation_file \
              --unique_id $library_id \
              --annotations_column "$annotations_column" \
              --celltype_value "$celltype_value" \
              --metaprograms_file $skc_metaprograms_file \
              --output_file $scores_output_file
        fi
    done
done
