#!/bin/bash
#SBATCH -J TriangleExperiment               # Job name
#SBATCH -N 1                   # Number of nodes
#SBATCH -n 1                  # Number of tasks
#SBATCH --time 24:00:00
#SBATCH --mem=8G
#SBATCH --constraint=broadwell
#SBATCH -o outputs/partition_%A_%a.txt       # Standard output file
#SBATCH -e errors/partition_%A_%a.txt        # Standard error file

module load gurobi
srun python partition_experiment.py $SLURM_ARRAY_TASK_ID 2 2
