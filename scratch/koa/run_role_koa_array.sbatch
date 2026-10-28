#!/bin/bash
#SBATCH --job-name=role_sbi_arr
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=1500M
#SBATCH --time=01:00:00
#SBATCH --array=1-20
#SBATCH --output=%x_%A_%a.out

# The 10,000 models as a Slurm array of small jobs instead of one 48-core job,
# because small jobs start much sooner. The chunk size is 10 models per core, so
# every size of job does the same amount of work per chunk (about 8 minutes at
# most), and each array task runs about 5 chunks (about 40 minutes at most).
# Number of array tasks for 10,000 models = 200 / cores per task:
#   10 cores: --array=1-20   4 cores: --array=1-50
#    2 cores: --array=1-100  1 core:  --array=1-200
#
# Submit from the folder that holds this file and run_role_koa.R:
#   sbatch run_role_koa_array.sbatch
# Other sizes, set on the command line (they override the lines above), e.g.
#   sbatch --cpus-per-task=4 --array=1-50 run_role_koa_array.sbatch
# Check estimated start times without submitting:
#   sbatch --test-only --cpus-per-task=4 --array=1-50 run_role_koa_array.sbatch
# Use --partition=shared,kill-shared to let Slurm pick whichever starts first
# (kill-shared jobs can be preempted; finished chunks are kept).
#
# Use the SAME cores-per-task for every submission that writes to the same
# output folder, because the chunk size depends on it.
# If tasks are killed or time out, submit the same command again; finished
# chunks are skipped. Results are combined into results_all.rds by whichever
# task finishes the last chunk.

# one thread per worker so the workers don't fight over cores
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

module load lang/R/4.5.2-gfbf-2025b
module load lang/Anaconda3/2024.02-1

# 10 models per core per chunk
CHUNK=$(( SLURM_CPUS_PER_TASK * 10 ))

cd "$SLURM_SUBMIT_DIR"
Rscript run_role_koa.R 10000 "$CHUNK" role_koa_out_10000_array
