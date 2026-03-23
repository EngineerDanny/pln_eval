#!/bin/bash
RSCRIPT=/projects/genomic-ml/da2343/soak-r/bin/Rscript
SCRIPT=/projects/genomic-ml/da2343/PLN/pln_eval/model_comp_subsample.R

# Subsample large datasets at N=100, 500, 1000, 2000
# Usage: model_comp_subsample.R <path> <n_samples> [memory_mb] [walltime_hours] [seed]

# deblur family (full: 9511 x 257)
srun --mem=2G $RSCRIPT $SCRIPT family/american_gut_1250/deblur_125nt_no_blooms_family.tsv.gz 100 2048 1 &
srun --mem=2G $RSCRIPT $SCRIPT family/american_gut_1250/deblur_125nt_no_blooms_family.tsv.gz 500 2048 1 &
srun --mem=4G $RSCRIPT $SCRIPT family/american_gut_1250/deblur_125nt_no_blooms_family.tsv.gz 1000 4096 2 &
srun --mem=4G $RSCRIPT $SCRIPT family/american_gut_1250/deblur_125nt_no_blooms_family.tsv.gz 2000 4096 2 &

# mbqc family (full: 18364 x 257)
srun --mem=2G $RSCRIPT $SCRIPT family/mbqc/mbqc_integrated_otus_family.tsv.gz 100 2048 1 &
srun --mem=2G $RSCRIPT $SCRIPT family/mbqc/mbqc_integrated_otus_family.tsv.gz 500 2048 1 &
srun --mem=4G $RSCRIPT $SCRIPT family/mbqc/mbqc_integrated_otus_family.tsv.gz 1000 4096 2 &
srun --mem=4G $RSCRIPT $SCRIPT family/mbqc/mbqc_integrated_otus_family.tsv.gz 2000 4096 2 &

# emp family (full: 27398 x 493)
srun --mem=2G $RSCRIPT $SCRIPT family/emp_release1_closed_ref_gg/emp_cr_gg_13_8.release1_family.tsv.gz 100 2048 1 &
srun --mem=2G $RSCRIPT $SCRIPT family/emp_release1_closed_ref_gg/emp_cr_gg_13_8.release1_family.tsv.gz 500 2048 1 &
srun --mem=4G $RSCRIPT $SCRIPT family/emp_release1_closed_ref_gg/emp_cr_gg_13_8.release1_family.tsv.gz 1000 4096 2 &
srun --mem=4G $RSCRIPT $SCRIPT family/emp_release1_closed_ref_gg/emp_cr_gg_13_8.release1_family.tsv.gz 2000 4096 2 &

wait
echo "All subsample submissions complete."
