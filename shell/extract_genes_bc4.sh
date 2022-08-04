#!/bin/bash

#SBATCH --array=1-4
#SBATCH --job-name=extract_genes
#SBATCH --partition=veryshort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --mem=32GB
#SBATCH --output=/user/work/tb13101/epi_gen_comp/slurm-output/extract_genes_%a.out
#SBATCH --error=/user/work/tb13101/epi_gen_comp/slurm-output/extract_genes_%a.error

scratch_dir="/user/work/tb13101/epi_gen_comp"
cd $scratch_dir
home_dir="/user/home/tb13101/projects/epi_gen_comp"
script=${home_dir}"/R/extract_genes_for_enrichment.R"
param_row=${SLURM_ARRAY_TASK_ID}

## Run analysis
srun Rscript $script ${param_row} "" $home_dir