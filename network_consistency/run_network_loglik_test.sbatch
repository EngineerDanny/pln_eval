#!/bin/bash
#SBATCH --job-name=network_loglik
#SBATCH --partition=core
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=/projects/genomic-ml/da2343/PLN/pln_eval/out/slurm/network_loglik_test_%A_%a.out
#SBATCH --error=/projects/genomic-ml/da2343/PLN/pln_eval/out/slurm/network_loglik_test_%A_%a.err

set -euo pipefail

BASE=/projects/genomic-ml/da2343/PLN/pln_eval
RSCRIPT=/projects/genomic-ml/da2343/soak-r/bin/Rscript
SCRIPT=$BASE/network_consistency/network_loglik_cv_worker.R

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "This wrapper requires an explicit array range."
  echo "Example: sbatch --array=1-12 $0 amgut2_update.csv"
  exit 2
fi

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch --array=1-N $0 <dataset-or-dataset-list>"
  exit 2
fi

DATASET=$1
TASK_ID=${SLURM_ARRAY_TASK_ID}

mkdir -p "$BASE/out/slurm"
cd "$BASE"

echo "Starting network_loglik_cv_worker.R at $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Dataset: $DATASET"
echo "Array task: $TASK_ID"

"$RSCRIPT" "$SCRIPT" "$DATASET" "$TASK_ID"

echo "Finished at $(date)"
