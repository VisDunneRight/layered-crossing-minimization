#!/bin/bash
#SBATCH -J Pictures               # Job name
#SBATCH -N 1                   # Number of nodes
#SBATCH -n 1                  # Number of tasks
#SBATCH --time=6:00:00
#SBATCH --mem=8G
#SBATCH --constraint=broadwell
#SBATCH -o outputs/pics_%A_%a.txt       # Standard output file
#SBATCH -e errors/pics_%A_%a.txt        # Standard error file

module load gurobi
srun python generate_pictures_all_metrics.py $SLURM_ARRAY_TASK_ID
