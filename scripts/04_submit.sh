#!/bin/bash
#SBATCH --job-name=grnboost2
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --qos=tc_normal_short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm_%x_%j.out
#SBATCH --error=logs/slurm_%x_%j.err

# Usage: sbatch scripts/04_submit.sh <celltype> [n_workers]
# Example (probe): sbatch scripts/04_submit.sh mature-endodermis 64

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

CELLTYPE="${1:?Usage: sbatch 04_submit.sh <celltype> [n_workers]}"
N_WORKERS="${2:-64}"

module load Miniconda3/25.11.1-1
source activate ~/.conda/envs/grn

echo "Job $SLURM_JOB_ID on $SLURM_NODELIST, celltype=$CELLTYPE, n_workers=$N_WORKERS"

python scripts/04_run_grnboost2.py \
  --celltype "$CELLTYPE" \
  --n-workers "$N_WORKERS" \
  --threads-per-worker 1
