#!/bin/bash
RSCRIPT=/projects/genomic-ml/da2343/soak-r/bin/Rscript
SCRIPT=/projects/genomic-ml/da2343/PLN/pln_eval/model_comp.R

# Small datasets (N<500): 2GB, 1h
srun --mem=2G $RSCRIPT $SCRIPT genus/cfmd/BertuzziAS_2018_genus.tsv.gz 2048 1 &
srun --mem=2G $RSCRIPT $SCRIPT family/microbiomehd/nash_ob_baker_family.tsv.gz 2048 1 &
srun --mem=2G $RSCRIPT $SCRIPT species/curated_metagenomic_data/20170526.HMP_2012.metaphlan_bugs_list.nasalcavity_species.tsv.gz 2048 1 &
srun --mem=2G $RSCRIPT $SCRIPT species/curated_metagenomic_data/HMP_2012.metaphlan_bugs_list.anterior_nares_species.tsv.gz 2048 1 &
srun --mem=2G $RSCRIPT $SCRIPT genus/curated_metagenomic_data/20181025.YuJ_2015.metaphlan_bugs_list.stool_genus.tsv.gz 2048 1 &
srun --mem=2G $RSCRIPT $SCRIPT genus/curated_metagenomic_data/20190424.CosteaPI_2017.metaphlan_bugs_list.stool_genus.tsv.gz 2048 1 &
srun --mem=2G $RSCRIPT $SCRIPT genus/microbiomehd/cdi_schubert_genus.tsv.gz 2048 1 &
srun --mem=2G $RSCRIPT $SCRIPT genus/curated_metagenomic_data/2021-03-31.VilaAV_2018.relative_abundance_genus.tsv.gz 2048 1 &

# Medium datasets (N>500): 4GB, 2h
srun --mem=4G $RSCRIPT $SCRIPT family/microbiomehd/ob_zeevi_family.tsv.gz 4096 2 &
srun --mem=4G $RSCRIPT $SCRIPT family/diabimmune_three_country/diabimmune_karelia_16s_family.tsv.gz 4096 2 &

wait
echo "All submissions complete."
