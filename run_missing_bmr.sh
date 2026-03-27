#!/bin/bash
set -euo pipefail

RSCRIPT=/projects/genomic-ml/da2343/soak-r/bin/Rscript
MODEL_SCRIPT=/projects/genomic-ml/da2343/PLN/pln_eval/model_comp.R
COLLECT_SCRIPT=/projects/genomic-ml/da2343/PLN/pln_eval/collect_wait.R
BASE_DIR=/projects/genomic-ml/da2343/PLN/pln_eval
SLURM_DIR=$BASE_DIR/out/slurm

mkdir -p "$SLURM_DIR"

prev_collect_job=

submit_dataset() {
  local dataset=$1
  local outer_mem=$2
  local outer_time=$3
  local inner_mem_mb=$4
  local inner_walltime_h=$5
  local safe_name model_job collect_job
  local -a dep_args=()

  safe_name=$(printf '%s' "$dataset" | tr -c 'A-Za-z0-9_.-' '_' | cut -c1-100)

  if [[ -n "$prev_collect_job" ]]; then
    dep_args=(--dependency "afterok:$prev_collect_job")
  fi

  model_job=$(
    sbatch --parsable "${dep_args[@]}" \
      --job-name "mc_$safe_name" \
      --cpus-per-task 1 \
      --mem "$outer_mem" \
      --time "$outer_time" \
      --output "$SLURM_DIR/model_comp_submit_%j.out" \
      --error "$SLURM_DIR/model_comp_submit_%j.err" \
      --wrap "cd $BASE_DIR && $RSCRIPT $MODEL_SCRIPT $dataset $inner_mem_mb $inner_walltime_h"
  )

  collect_job=$(
    sbatch --parsable \
      --dependency "afterok:$model_job" \
      --job-name "collect_$safe_name" \
      --cpus-per-task 1 \
      --mem 4G \
      --time 24:00:00 \
      --output "$SLURM_DIR/collect_wait_%j.out" \
      --error "$SLURM_DIR/collect_wait_%j.err" \
      --wrap "cd $BASE_DIR && $RSCRIPT $COLLECT_SCRIPT $dataset"
  )

  prev_collect_job=$collect_job
  echo "$dataset model_job=$model_job collect_job=$collect_job"
}

submit_small() {
  submit_dataset "$1" 2G 01:00:00 2048 1
}

submit_medium() {
  submit_dataset "$1" 4G 01:30:00 4096 2
}

# Small-first pilot set: diverse body sites, diseases, taxonomic levels, and one CSV-backed dataset.
submit_small amgut1
submit_small 20170526.AsnicarF_2017.metaphlan_bugs_list.milk_genus
submit_small 20170526.Castro-NallarE_2015.metaphlan_bugs_list.oralcavity_genus
submit_small 20170526.ChngKR_2016.metaphlan_bugs_list.skin_genus
submit_small 20181025.OlmMR_2017.metaphlan_bugs_list.stool_genus
submit_small cdi_schubert_family
submit_small t1d_alkanani_family
submit_small ibd_gevers_2014_family
submit_small ob_ross_genus
submit_small nash_chan_family
submit_small crc_zhao_genus
submit_small hiv_lozupone_family
submit_small par_scheperjans_family
submit_medium 2021-10-14.MehtaRS_2018.relative_abundance_genus
submit_medium 20170526.AsnicarF_2017.metaphlan_bugs_list.milk_species

echo "Submitted 15 datasets in a dependency chain."
echo "Last collect job: $prev_collect_job"
