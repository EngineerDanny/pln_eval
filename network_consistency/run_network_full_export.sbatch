#!/bin/bash
#SBATCH --job-name=network_full_export
#SBATCH --partition=core
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=/projects/genomic-ml/da2343/PLN/pln_eval/out/slurm/network_full_export_%j.out
#SBATCH --error=/projects/genomic-ml/da2343/PLN/pln_eval/out/slurm/network_full_export_%j.err

set -euo pipefail

BASE=/projects/genomic-ml/da2343/PLN/pln_eval
RSCRIPT=/projects/genomic-ml/da2343/soak-r/bin/Rscript
SCRIPT=$BASE/network_consistency/network_fit_full_export.R

DATASET=${1:-amgut2_update.csv}
TOP_N=${2:-all}
OUT_TAG=${3:-}

mkdir -p "$BASE/out/slurm"
cd "$BASE"

echo "Starting network_fit_full_export.R at $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Dataset: $DATASET"
echo "Top taxa: $TOP_N"
if [[ -n "$OUT_TAG" ]]; then
  echo "Output tag: $OUT_TAG"
  "$RSCRIPT" "$SCRIPT" "$DATASET" "$TOP_N" "$OUT_TAG"
else
  "$RSCRIPT" "$SCRIPT" "$DATASET" "$TOP_N"
fi

echo "Finished at $(date)"
