#!/bin/bash
#SBATCH --job-name=clr_sequence
#SBATCH --time=48:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --chdir=/projects/genomic-ml/da2343/PLN/pln_eval
#SBATCH --output=/projects/genomic-ml/da2343/PLN/pln_eval/logs/clr_sequence_%j.out
#SBATCH --error=/projects/genomic-ml/da2343/PLN/pln_eval/logs/clr_sequence_%j.err

set -euo pipefail

BASE_DIR=/projects/genomic-ml/da2343/PLN/pln_eval
RSCRIPT=/projects/genomic-ml/da2343/soak-r/bin/Rscript
MODEL=$BASE_DIR/clr_sensitivity/model_comp_clr.R
COLLECT=$BASE_DIR/clr_sensitivity/collect_results_clr.R

run() {
  local dataset=$1
  local memory_mb=$2
  local walltime_hours=$3

  echo "Submitting CLR sensitivity analysis for $dataset"
  "$RSCRIPT" "$MODEL" "$dataset" "$memory_mb" "$walltime_hours"

  echo "Waiting for $dataset to finish"
  "$RSCRIPT" "$COLLECT" "$dataset"
  echo "Completed $dataset"
}

run amgut2 4096 2
run crc_zeller 4096 2
run hiv_lozupone_family 2048 1
run cdi_schubert_family 2048 1
run diabimmune_karelia_16s_family 4096 2
run amgut1 2048 1
run 20190424.CosteaPI_2017.metaphlan_bugs_list.stool_genus 2048 1
run mbqc_integrated_otus_family 8192 4
