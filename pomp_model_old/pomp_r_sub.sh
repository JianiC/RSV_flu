#!/bin/bash

#SBATCH --job-name=new_null
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=3gb
#SBATCH --cpus-per-task=10
#SBATCH --time=150:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err




cd $SLURM_SUBMIT_DIR


module load Miniconda3/4.7.10

source activate /home/jc54391/conda_R_env2

R CMD BATCH new_pompmodel_null.R
