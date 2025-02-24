# Benchmarking Instructions
## Running the Benchmark in Series
To run the benchmark on your local machine:
1. Navigate to the design_constraints directory in your terminal
2. Run `python3 benchmark_exec.py`

Or, open `benchmark_exec.py` in your python IDE of choice and use the run feature

## Running the Benchmark in Parallel
To run the benchmark in parallel on your local machine:
1. Navigate to the design_constraints directory in your terminal
2. Run `bash benchmark_parallel.sh`
## Running the Benchmark Remotely using SLURM
This is recommended for large benchmark experiments with many combinations of conditions. To run the benchmark in parallel on a remote cluster that uses SLURM:
1. Access the cluster and navigate to the design_constraints directory
2. If necessary, activate a python virtual environment using as conda or pip with all required packages installed
3. If necessary, modify `benchmark_slurm.sh` to load any needed modules prior to the `srun` call, e.g. `module load gurobi`4. Run the following: `sbatch --array=0-11 benchmark_slurm.sh`

*If you need to run the benchmark in series using SLURM, use the steps above with the following modifications:
- delete $SLURM_ARRAY_TASK_ID from benchmark_slurm.sh
- run `sbatch benchmark_slurm.sh`
