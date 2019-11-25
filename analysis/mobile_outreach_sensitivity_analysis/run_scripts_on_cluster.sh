#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --time=2-1:00:00
#SBATCH --mem-per-cpu=15GB
#SBATCH --job-name=simulate_interventions
#SBATCH --error=sbatch-out/%j.err
#SBATCH --output=sbatch-out/%j.out

source ~/load-gcc-and-R.sh
srun Rscript simulate_interventions.R $1
echo $1
