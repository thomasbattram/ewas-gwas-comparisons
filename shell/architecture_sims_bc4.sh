#!/bin/bash

#SBATCH --array=1-100
#SBATCH --job-name=arch_sims
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=70:00:00
#SBATCH --mem=32GB
#SBATCH --output=/user/work/tb13101/epi_gen_comp/slurm-output/arch_sims_%a.out
#SBATCH --error=/user/work/tb13101/epi_gen_comp/slurm-output/arch_sims_%a.error

scratch_dir="/user/work/tb13101/epi_gen_comp"
cd $scratch_dir
home_dir="/user/home/tb13101/projects/epi_gen_comp"
script=${home_dir}"/R/architecture_sims.R"
param_row=${SLURM_ARRAY_TASK_ID}

## Run analysis
srun Rscript $script ${param_row} "egfr" $home_dir
srun Rscript $script ${param_row} "urate" $home_dir