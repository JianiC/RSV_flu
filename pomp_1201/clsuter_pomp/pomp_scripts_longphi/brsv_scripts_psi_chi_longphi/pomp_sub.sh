#!/bin/bash

#SBATCH --job-name=5_brsv_chi
#SBATCH --partition=bahl_salv_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --cpus-per-task=60 
#SBATCH --time=150:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --constraint=EDR




cd $SLURM_SUBMIT_DIR


module load R/4.0.0-foss-2019b

R CMD BATCH ./hhs5_brsv_chi_script.R