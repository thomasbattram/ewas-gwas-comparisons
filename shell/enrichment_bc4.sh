#!/bin/bash

#SBATCH --array=14-15
#SBATCH --job-name=enrich_tests
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=30:00:00
#SBATCH --mem=32GB
#SBATCH --output=/user/work/tb13101/epi_gen_comp/slurm-output/enrich_tests_%a.out
#SBATCH --error=/user/work/tb13101/epi_gen_comp/slurm-output/enrich_tests_%a.error

scratch_dir="/user/work/tb13101/epi_gen_comp"
cd $scratch_dir
home_dir="/user/home/tb13101/projects/epi_gen_comp"
script=${home_dir}"/R/enrichment_tests.R"
param_row=${SLURM_ARRAY_TASK_ID}

## Run analysis
srun Rscript $script ${param_row} "" $home_dir