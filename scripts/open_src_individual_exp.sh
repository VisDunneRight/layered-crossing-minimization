#!/bin/bash
#SBATCH -J OpenSrcIndExperiment               # Job name
#SBATCH -N 1                   # Number of nodes
#SBATCH -n 1                  # Number of tasks
#SBATCH --time 24:00:00
#SBATCH --constraint=broadwell
#SBATCH --mem=8G
#SBATCH -o outputs/opensrc_%A_%a.txt       # Standard output file
#SBATCH -e errors/opensrc_%A_%a.txt        # Standard error file

srun python experiments.py 5 $SLURM_ARRAY_TASK_ID