#!/bin/bash
#SBATCH -J OpenSrcIndExperiment               # Job name
#SBATCH -N 1                   # Number of nodes
#SBATCH -n 1                  # Number of tasks
#SBATCH --time 12:00:00
#SBATCH -o output_%A_%a.txt       # Standard output file
#SBATCH -e error_%A_%a.txt        # Standard error file

srun python ../experiments.py 5 $SLURM_ARRAY_TASK_ID
