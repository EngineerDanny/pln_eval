## Compute mean absolute correlation (MAC) for all full-sample benchmarked datasets

data_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval/data"

load_counts_matrix <- function(tsv_path) {
  raw <- read.table(gzfile(tsv_path), sep = "\t", header = TRUE,
                    row.names = 1, check.names = FALSE)
  counts <- t(as.matrix(raw))
  counts <- counts[rowSums(counts) > 0, colSums(counts) > 0]
  counts
}

compute_mac <- function(counts) {
  N <- nrow(counts)
  p <- ncol(counts)
  X <- log1p(counts)
  cor_mat <- cor(X)
  diag(cor_mat) <- 0
  mac <- mean(abs(cor_mat), na.rm = TRUE)
  data.frame(N = N, p = p, N_over_p = round(N / p, 2),
             mac = round(mac, 3), stringsAsFactors = FALSE)
}

# Full-sample datasets only
datasets <- list(
  list(name = "AsnicarF_2017 species",
       path = "species/curated_metagenomic_data/2021-10-14.AsnicarF_2017.relative_abundance_species.tsv.gz"),
  list(name = "HMP_2012 nasalcavity species",
       path = "species/curated_metagenomic_data/20170526.HMP_2012.metaphlan_bugs_list.nasalcavity_species.tsv.gz"),
  list(name = "HMP_2012 anterior_nares species",
       path = "species/curated_metagenomic_data/HMP_2012.metaphlan_bugs_list.anterior_nares_species.tsv.gz"),
  list(name = "PetersBA_2019 genus",
       path = "genus/curated_metagenomic_data/2022-04-13.PetersBA_2019.relative_abundance_genus.tsv.gz"),
  list(name = "OlmMR_2017 genus",
       path = "genus/curated_metagenomic_data/2021-03-31.OlmMR_2017.relative_abundance_genus.tsv.gz"),
  list(name = "KarlssonFH_2013 stool genus",
       path = "genus/curated_metagenomic_data/20180425.KarlssonFH_2013.metaphlan_bugs_list.stool_genus.tsv.gz"),
  list(name = "YuJ_2015 stool genus",
       path = "genus/curated_metagenomic_data/20181025.YuJ_2015.metaphlan_bugs_list.stool_genus.tsv.gz"),
  list(name = "CosteaPI_2017 stool genus",
       path = "genus/curated_metagenomic_data/20190424.CosteaPI_2017.metaphlan_bugs_list.stool_genus.tsv.gz"),
  list(name = "MehtaRS_2018 genus",
       path = "genus/curated_metagenomic_data/2021-03-31.MehtaRS_2018.relative_abundance_genus.tsv.gz"),
  list(name = "ShaoY_2019 genus",
       path = "genus/curated_metagenomic_data/2021-03-31.ShaoY_2019.relative_abundance_genus.tsv.gz"),
  list(name = "VilaAV_2018 genus",
       path = "genus/curated_metagenomic_data/2021-03-31.VilaAV_2018.relative_abundance_genus.tsv.gz"),
  list(name = "RubelMA_2020 species",
       path = "species/curated_metagenomic_data/2021-10-14.RubelMA_2020.relative_abundance_species.tsv.gz"),
  list(name = "ra_littman family",
       path = "family/microbiomehd/ra_littman_family.tsv.gz"),
  list(name = "otu_table_psn_v13 family",
       path = "family/hmp16s/otu_table_psn_v13_family.tsv.gz"),
  list(name = "ob_zeevi family",
       path = "family/microbiomehd/ob_zeevi_family.tsv.gz"),
  list(name = "nash_ob_baker family",
       path = "family/microbiomehd/nash_ob_baker_family.tsv.gz"),
  list(name = "diabimmune_karelia family",
       path = "family/diabimmune_three_country/diabimmune_karelia_16s_family.tsv.gz"),
  list(name = "crc_baxter",
       path = "family/microbiomehd/crc_baxter_family.tsv.gz"),
  list(name = "crc_zeller",
       path = "family/microbiomehd/crc_zeller_family.tsv.gz"),
  list(name = "cdi_schubert genus",
       path = "genus/microbiomehd/cdi_schubert_genus.tsv.gz"),
  list(name = "BertuzziAS_2018 genus",
       path = "genus/cfmd/BertuzziAS_2018_genus.tsv.gz")
)

results <- list()
for (ds in datasets) {
  tsv_path <- file.path(data_dir, ds$path)
  if (!file.exists(tsv_path)) {
    cat("SKIP (not found):", ds$name, "-", tsv_path, "\n")
    next
  }
  cat("Processing:", ds$name, "...")
  counts <- load_counts_matrix(tsv_path)
  metrics <- compute_mac(counts)
  metrics$dataset <- ds$name
  results[[length(results) + 1]] <- metrics
  cat(" done (N=", metrics$N, "p=", metrics$p, "mac=", metrics$mac, ")\n")
}

result_df <- do.call(rbind, results)
result_df <- result_df[, c("dataset", "N", "p", "N_over_p", "mac")]

cat("\n\n=== MEAN ABSOLUTE CORRELATION ===\n\n")
fmt <- "%-35s %5s %5s %6s %5s\n"
cat(sprintf(fmt, "Dataset", "N", "p", "N/p", "MAC"))
cat(paste(rep("-", 60), collapse = ""), "\n")
for (i in 1:nrow(result_df)) {
  r <- result_df[i, ]
  cat(sprintf(fmt, r$dataset, r$N, r$p, r$N_over_p, r$mac))
}

out_path <- "/projects/genomic-ml/da2343/PLN/pln_eval/covariance_metrics.csv"
write.csv(result_df, out_path, row.names = FALSE)
cat("\nSaved to:", out_path, "\n")
