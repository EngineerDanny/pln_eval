#!/bin/bash
#SBATCH --job-name=network_loglik_collect
#SBATCH --partition=core
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=/projects/genomic-ml/da2343/PLN/pln_eval/out/slurm/network_loglik_collect_%j.out
#SBATCH --error=/projects/genomic-ml/da2343/PLN/pln_eval/out/slurm/network_loglik_collect_%j.err

set -euo pipefail

BASE=/projects/genomic-ml/da2343/PLN/pln_eval
RSCRIPT=/projects/genomic-ml/da2343/soak-r/bin/Rscript
SCRIPT=$BASE/network_consistency/network_loglik_cv_collect.R

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch $0 <dataset-or-dataset-list>"
  exit 2
fi

DATASET=$1

mkdir -p "$BASE/out/slurm"
cd "$BASE"

echo "Starting network_loglik_cv_collect.R at $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Dataset: $DATASET"

"$RSCRIPT" "$SCRIPT" "$DATASET"

echo "Finished at $(date)"
