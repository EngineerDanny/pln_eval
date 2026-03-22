#!/bin/bash
# Submit model_comp.R for diverse dataset subset
# Usage: Rscript model_comp.R <path_or_name> [memory_mb] [walltime_hours]

RSCRIPT=/projects/genomic-ml/da2343/soak-r/bin/Rscript
SCRIPT=/projects/genomic-ml/da2343/PLN/pln_eval/model_comp.R

# Small datasets (N<200): 2GB, 1h
srun --mem=2G $RSCRIPT $SCRIPT species/curated_metagenomic_data/2021-10-14.AsnicarF_2017.relative_abundance_species.tsv.gz 2048 1 &
srun --mem=2G $RSCRIPT $SCRIPT genus/curated_metagenomic_data/2022-04-13.PetersBA_2019.relative_abundance_genus.tsv.gz 2048 1 &
srun --mem=2G $RSCRIPT $SCRIPT genus/curated_metagenomic_data/2021-03-31.OlmMR_2017.relative_abundance_genus.tsv.gz 2048 1 &
srun --mem=2G $RSCRIPT $SCRIPT family/microbiomehd/ra_littman_family.tsv.gz 2048 1 &
srun --mem=2G $RSCRIPT $SCRIPT genus/curated_metagenomic_data/20180425.KarlssonFH_2013.metaphlan_bugs_list.stool_genus.tsv.gz 2048 1 &
srun --mem=2G $RSCRIPT $SCRIPT species/curated_metagenomic_data/2021-10-14.RubelMA_2020.relative_abundance_species.tsv.gz 2048 1 &

# Medium datasets (N=928-1644): 4GB, 2h
srun --mem=4G $RSCRIPT $SCRIPT genus/curated_metagenomic_data/2021-03-31.MehtaRS_2018.relative_abundance_genus.tsv.gz 4096 2 &
srun --mem=4G $RSCRIPT $SCRIPT genus/curated_metagenomic_data/2021-03-31.ShaoY_2019.relative_abundance_genus.tsv.gz 4096 2 &

# Large dataset (N=2910): 4GB, 4h
srun --mem=4G $RSCRIPT $SCRIPT family/hmp16s/otu_table_psn_v13_family.tsv.gz 4096 4 &

wait
echo "All submissions complete."
