#!/bin/bash
#SBATCH -J FixedDensityExperiment               # Job name
#SBATCH -N 1                   # Number of nodes
#SBATCH -n 1                  # Number of tasks
#SBATCH -o outputs/output_%A_%a.txt       # Standard output file
#SBATCH -e errors/error_%A_%a.txt        # Standard error file

module load gurobi
srun python experiments.py 4 $SLURM_ARRAY_TASK_ID
