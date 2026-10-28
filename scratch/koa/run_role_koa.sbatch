#!/bin/bash
#SBATCH --job-name=role_sbi
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=1500M
#SBATCH --time=06:00:00
#SBATCH --output=%x_%j.out

# Submit from the folder that holds both this file and run_role_koa.R:
#   sbatch run_role_koa.sbatch
# Do a small test first with run_role_koa_sandbox.sbatch (16 models).
# Sizing comes from that test (job 15158122: 8 workers, 2.9 GB peak, about
# 1 minute per chunk of 8 models); the full run is 21 chunks of 480 models.
# If the job is killed, submit the same command again to resume.

# one thread per worker so the workers don't fight over cores
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# load the same R and conda that roleR and the r-roleR env were installed with;
# roleR finds the conda env (in ~/.conda/envs) by itself
module load lang/R/4.5.2-gfbf-2025b
module load lang/Anaconda3/2024.02-1

cd "$SLURM_SUBMIT_DIR"
Rscript run_role_koa.R "${N_MODELS:-10000}" "${CHUNK_SIZE:-480}"
