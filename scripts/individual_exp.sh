#!/bin/bash
#SBATCH -J IndividualSwitchExperiment               # Job name
#SBATCH -N 1                   # Number of nodes
#SBATCH -n 1                  # Number of tasks
#SBATCH -o output_%A_%a.txt       # Standard output file
#SBATCH -e error_%A_%a.txt        # Standard error file

module load gurobi
srun python ../experiments.py 1 $SLURM_ARRAY_TASK_ID
