#!/bin/bash

#SBATCH --job-name=jobname
#SBATCH --partition=bahl_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80gb
#SBATCH --cpus-per-task=93 
#SBATCH --time=150:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --constraint=EDR




cd $SLURM_SUBMIT_DIR


module load R/4.0.0-foss-2019b

