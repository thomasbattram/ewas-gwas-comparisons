#!/bin/bash

#SBATCH --array=1-100
#SBATCH --job-name=method_sims
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=16GB
#SBATCH --output=/user/work/tb13101/epi_gen_comp/slurm-output/method_sims_%a.out
#SBATCH --error=/user/work/tb13101/epi_gen_comp/slurm-output/method_sims_%a.error

scratch_dir="/user/work/tb13101/epi_gen_comp"
cd $scratch_dir
home_dir="/user/home/tb13101/projects/epi_gen_comp"
script=${home_dir}"/R/comparing_overlap_methods_sims_gwas_power.R"
split=${SLURM_ARRAY_TASK_ID}
sim_n="2"
genes_or_proteins="genes"
database_arg="go"

## Run analysis
srun Rscript $script ${sim_n} ${split} ${genes_or_proteins} ${database_arg} ${home_dir}
