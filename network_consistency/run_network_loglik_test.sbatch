#!/bin/bash
#SBATCH --job-name=network_loglik
#SBATCH --partition=core
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=/projects/genomic-ml/da2343/PLN/pln_eval/out/slurm/network_loglik_test_%j.out
#SBATCH --error=/projects/genomic-ml/da2343/PLN/pln_eval/out/slurm/network_loglik_test_%j.err

set -euo pipefail

BASE=/projects/genomic-ml/da2343/PLN/pln_eval
RSCRIPT=/projects/genomic-ml/da2343/soak-r/bin/Rscript
SCRIPT=$BASE/network_consistency/network_loglik_test.R
DATASET=${1:-amgut2_update.csv}

mkdir -p "$BASE/out/slurm"
cd "$BASE"

echo "Starting network_loglik_test.R at $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Dataset: $DATASET"

"$RSCRIPT" "$SCRIPT" "$DATASET"

echo "Finished at $(date)"
