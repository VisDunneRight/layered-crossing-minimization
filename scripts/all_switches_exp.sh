#!/bin/bash
#SBATCH -J AllSwitchesExperiment               # Job name
#SBATCH -N 1                   # Number of nodes
#SBATCH -n 1                  # Number of tasks
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --constraint=broadwell
#SBATCH -o outputs/all_%A_%a.txt       # Standard output file
#SBATCH -e errors/all_%A_%a.txt        # Standard error file

module load gurobi
srun python experiments.py 2 $SLURM_ARRAY_TASK_ID
