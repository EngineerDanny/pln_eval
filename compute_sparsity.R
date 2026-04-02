## Compute sparsity for all 34 benchmarked datasets

data_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval/data"
out_dir  <- "/projects/genomic-ml/da2343/PLN/pln_eval/out"

load_counts_tsv <- function(path) {
  raw    <- read.table(gzfile(path), sep = "\t", header = TRUE,
                       row.names = 1, check.names = FALSE)
  counts <- t(as.matrix(raw))
  counts <- counts[rowSums(counts) > 0, colSums(counts) > 0, drop = FALSE]
  counts
}

load_counts_csv <- function(path) {
  raw    <- read.csv(path, check.names = FALSE)
  counts <- as.matrix(raw)
  storage.mode(counts) <- "numeric"
  counts <- counts[rowSums(counts) > 0, colSums(counts) > 0, drop = FALSE]
  counts
}

sparsity <- function(counts) round(mean(counts == 0) * 100, 2)

# All dataset paths keyed by bmr CSV name
dataset_paths <- list(
  "2021-10-14.AsnicarF_2017.relative_abundance_species" =
    file.path(data_dir, "species/curated_metagenomic_data/2021-10-14.AsnicarF_2017.relative_abundance_species.tsv.gz"),
  "2022-04-13.PetersBA_2019.relative_abundance_genus" =
    file.path(data_dir, "genus/curated_metagenomic_data/2022-04-13.PetersBA_2019.relative_abundance_genus.tsv.gz"),
  "20170526.Castro-NallarE_2015.metaphlan_bugs_list.oralcavity_genus" =
    file.path(data_dir, "genus/curated_metagenomic_data/20170526.Castro-NallarE_2015.metaphlan_bugs_list.oralcavity_genus.tsv.gz"),
  "20170526.ChngKR_2016.metaphlan_bugs_list.skin_genus" =
    file.path(data_dir, "genus/curated_metagenomic_data/20170526.ChngKR_2016.metaphlan_bugs_list.skin_genus.tsv.gz"),
  "20170526.HMP_2012.metaphlan_bugs_list.nasalcavity_species" =
    file.path(data_dir, "species/curated_metagenomic_data/20170526.HMP_2012.metaphlan_bugs_list.nasalcavity_species.tsv.gz"),
  "HMP_2012.metaphlan_bugs_list.anterior_nares_species" =
    file.path(data_dir, "species/curated_metagenomic_data/HMP_2012.metaphlan_bugs_list.anterior_nares_species.tsv.gz"),
  "2021-10-14.RubelMA_2020.relative_abundance_species" =
    file.path(data_dir, "species/curated_metagenomic_data/2021-10-14.RubelMA_2020.relative_abundance_species.tsv.gz"),
  "crc_zeller" =
    file.path(data_dir, "family/microbiomehd/crc_zeller_family.tsv.gz"),
  "ob_ross_genus" =
    file.path(data_dir, "genus/microbiomehd/ob_ross_genus.tsv.gz"),
  "20181025.YuJ_2015.metaphlan_bugs_list.stool_genus" =
    file.path(data_dir, "genus/curated_metagenomic_data/20181025.YuJ_2015.metaphlan_bugs_list.stool_genus.tsv.gz"),
  "hiv_lozupone_family" =
    file.path(data_dir, "family/microbiomehd/hiv_lozupone_family.tsv.gz"),
  "20180425.KarlssonFH_2013.metaphlan_bugs_list.stool_genus" =
    file.path(data_dir, "genus/curated_metagenomic_data/20180425.KarlssonFH_2013.metaphlan_bugs_list.stool_genus.tsv.gz"),
  "nash_ob_baker_family" =
    file.path(data_dir, "family/microbiomehd/nash_ob_baker_family.tsv.gz"),
  "crc_zhao_genus" =
    file.path(data_dir, "genus/microbiomehd/crc_zhao_genus.tsv.gz"),
  "BertuzziAS_2018_genus" =
    file.path(data_dir, "genus/cfmd/BertuzziAS_2018_genus.tsv.gz"),
  "20190424.CosteaPI_2017.metaphlan_bugs_list.stool_genus" =
    file.path(data_dir, "genus/curated_metagenomic_data/20190424.CosteaPI_2017.metaphlan_bugs_list.stool_genus.tsv.gz"),
  "cdi_schubert_genus" =
    file.path(data_dir, "genus/microbiomehd/cdi_schubert_genus.tsv.gz"),
  "2021-03-31.VilaAV_2018.relative_abundance_genus" =
    file.path(data_dir, "genus/curated_metagenomic_data/2021-03-31.VilaAV_2018.relative_abundance_genus.tsv.gz"),
  "amgut2" =
    file.path(data_dir, "amgut2_update.csv"),
  "2021-03-31.OlmMR_2017.relative_abundance_genus" =
    file.path(data_dir, "genus/curated_metagenomic_data/2021-03-31.OlmMR_2017.relative_abundance_genus.tsv.gz"),
  "nash_chan_family" =
    file.path(data_dir, "family/microbiomehd/nash_chan_family.tsv.gz"),
  "amgut1" =
    file.path(data_dir, "amgut1_update.csv"),
  "ra_littman_family" =
    file.path(data_dir, "family/microbiomehd/ra_littman_family.tsv.gz"),
  "par_scheperjans_family" =
    file.path(data_dir, "family/microbiomehd/par_scheperjans_family.tsv.gz"),
  "ibd_gevers_2014_family" =
    file.path(data_dir, "family/microbiomehd/ibd_gevers_2014_family.tsv.gz"),
  "t1d_alkanani_family" =
    file.path(data_dir, "family/microbiomehd/t1d_alkanani_family.tsv.gz"),
  "cdi_schubert_family" =
    file.path(data_dir, "family/microbiomehd/cdi_schubert_family.tsv.gz"),
  "2021-03-31.MehtaRS_2018.relative_abundance_genus" =
    file.path(data_dir, "genus/curated_metagenomic_data/2021-03-31.MehtaRS_2018.relative_abundance_genus.tsv.gz"),
  "crc_baxter" =
    file.path(data_dir, "family/microbiomehd/crc_baxter_family.tsv.gz"),
  "2021-03-31.ShaoY_2019.relative_abundance_genus" =
    file.path(data_dir, "genus/curated_metagenomic_data/2021-03-31.ShaoY_2019.relative_abundance_genus.tsv.gz"),
  "ob_zeevi_family" =
    file.path(data_dir, "family/microbiomehd/ob_zeevi_family.tsv.gz"),
  "otu_table_psn_v13_family" =
    file.path(data_dir, "family/hmp16s/otu_table_psn_v13_family.tsv.gz"),
  "diabimmune_karelia_16s_family" =
    file.path(data_dir, "family/diabimmune_three_country/diabimmune_karelia_16s_family.tsv.gz"),
  "mbqc_integrated_otus_family" =
    file.path(data_dir, "family/mbqc/mbqc_integrated_otus_family.tsv.gz")
)

# Load comparison table for winner info
cmp <- read.csv(file.path(out_dir, "full_comparison_table_v2.csv"), stringsAsFactors = FALSE)

results <- list()
for (ds in names(dataset_paths)) {
  path <- dataset_paths[[ds]]
  cat("Processing:", ds, "...")
  if (!file.exists(path)) {
    cat(" SKIP (not found)\n")
    next
  }
  counts <- if (grepl("\\.csv$", path)) load_counts_csv(path) else load_counts_tsv(path)
  sp <- sparsity(counts)
  row <- cmp[cmp$dataset == ds, ]
  results[[length(results) + 1]] <- data.frame(
    dataset  = ds,
    N        = nrow(counts),
    p        = ncol(counts),
    N_over_p = round(nrow(counts) / ncol(counts), 2),
    sparsity = sp,
    winner   = if (nrow(row) > 0) row$winner else NA,
    pln_vs_lasso_pct = if (nrow(row) > 0) row$pln_vs_lasso_pct else NA,
    stringsAsFactors = FALSE
  )
  cat(" sparsity =", sp, "%\n")
}

out <- do.call(rbind, results)
out <- out[order(out$N_over_p), ]

cat("\n=== SPARSITY SUMMARY ===\n")
fmt <- "%-55s %5s %5s %6s %8s %6s %s\n"
cat(sprintf(fmt, "Dataset", "N", "p", "N/p", "Sparsity", "Winner", "PLN vs Lasso"))
cat(paste(rep("-", 110), collapse=""), "\n")
for (i in seq_len(nrow(out))) {
  r <- out[i, ]
  cat(sprintf(fmt, substr(r$dataset, 1, 55), r$N, r$p, r$N_over_p,
              paste0(r$sparsity, "%"), r$winner, r$pln_vs_lasso_pct))
}

write.csv(out, file.path(out_dir, "sparsity_summary.csv"), row.names = FALSE)
cat("\nSaved to:", file.path(out_dir, "sparsity_summary.csv"), "\n")
