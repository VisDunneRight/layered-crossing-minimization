#!/bin/bash
#SBATCH -J HeuristicExperiment               # Job name
#SBATCH -N 1                   # Number of nodes
#SBATCH -n 1                  # Number of tasks
#SBATCH --time 24:00:00
#SBATCH --mem=8G
#SBATCH --constraint=broadwell
#SBATCH -o outputs/heuristic_%A_%a.txt       # Standard output file
#SBATCH -e errors/heuristic_%A_%a.txt        # Standard error file

srun python heuristic_exp.py $SLURM_ARRAY_TASK_ID
