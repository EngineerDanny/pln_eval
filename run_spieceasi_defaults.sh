#!/bin/bash
# Submit one resumable SPIEC-EASI default-settings job per network dataset.

set -euo pipefail

BASE_DIR=/projects/genomic-ml/da2343/PLN/pln_eval
RSCRIPT=/projects/genomic-ml/da2343/soak-r/bin/Rscript
RESULT_DIR=$BASE_DIR/out/spieceasi_default_by_dataset
LOG_DIR=$BASE_DIR/logs/spieceasi_defaults
mkdir -p "$RESULT_DIR" "$LOG_DIR"

datasets=(
  omm12
  omm12_keystone_2023
  pairinterax
  butyrate_assembly_2021
  host_fitness_2018
)

for dataset in "${datasets[@]}"; do
  if [[ "$dataset" == "pairinterax" ]]; then
    walltime=12:00:00
    memory=8G
  else
    walltime=02:00:00
    memory=4G
  fi

  output=$RESULT_DIR/${dataset}.csv
  job_id=$(sbatch --parsable \
    --job-name="spiec_${dataset}" \
    --time="$walltime" \
    --mem="$memory" \
    --cpus-per-task=1 \
    --chdir="$BASE_DIR" \
    --output="$LOG_DIR/${dataset}_%j.out" \
    --error="$LOG_DIR/${dataset}_%j.err" \
    --wrap="PATH=/projects/genomic-ml/da2343/soak-r/bin:\$PATH SPIECEASI_OUTPUT=$output $RSCRIPT $BASE_DIR/spieceasi_f1.R $dataset")
  printf '%s\t%s\n' "$dataset" "$job_id"
done
