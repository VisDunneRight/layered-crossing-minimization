#!/bin/bash
#SBATCH -J OSAllExperiment               # Job name
#SBATCH -N 1                   # Number of nodes
#SBATCH -n 1                  # Number of tasks
#SBATCH --time 12:00:00
#SBATCH --constraint=broadwell
#SBATCH --mem=8G
#SBATCH -o outputs/os_all_%A_%a.txt       # Standard output file
#SBATCH -e errors/os_all_%A_%a.txt        # Standard error file

srun python experiments.py 6 $SLURM_ARRAY_TASK_ID
