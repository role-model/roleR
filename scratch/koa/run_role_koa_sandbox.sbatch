#!/bin/bash
#SBATCH --job-name=role_sbi_test
#SBATCH --partition=sandbox
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1500M
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out

# Small test run in the sandbox partition (4 hour limit): 16 models in 2 chunks
# of 8, on 8 cores, so each worker runs one model per chunk and the time per
# chunk is about the time of the slowest model. Use it to see how long models take and how much memory
# they use before submitting run_role_koa.sbatch for the full 10,000.
#
# Submit from the folder that holds this file and run_role_koa.R:
#   sbatch run_role_koa_sandbox.sbatch
#
# Afterwards, check what was actually used (MaxRSS is per task, in KB):
#   sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,ReqMem,AllocCPUS,State
# Output goes to role_koa_out_300/ and the log shows the time per chunk.

# one thread per worker so the workers don't fight over cores
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# load the same R and conda that roleR and the r-roleR env were installed with;
# roleR finds the conda env (in ~/.conda/envs) by itself
module load lang/R/4.5.2-gfbf-2025b
module load lang/Anaconda3/2024.02-1

cd "$SLURM_SUBMIT_DIR"
Rscript run_role_koa.R 300 16
