#!/bin/bash
#SBATCH -J aesthetics               # Job name
#SBATCH --partition=short        # Use short partition (24hrs max)
#SBATCH -N 1                   # Number of nodes
#SBATCH -n 1                  # Number of tasks
#SBATCH --time 12:00:00           # Request 12 hours of compute time
#SBATCH --mem=8G           # Request 8GB memory
#SBATCH --constraint=broadwell           # Run job on the broadwell hardware cluster
#SBATCH -o outputs/aesthetics_%A_%a.txt       # Standard output file
#SBATCH -e errors/aesthetics_%A_%a.txt        # Standard error file

srun python aesthetics_exp.py $SLURM_ARRAY_TASK_ID
