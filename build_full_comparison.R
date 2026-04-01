## Build comprehensive comparison table from all bmr/ benchmark results
## Outputs: /projects/genomic-ml/da2343/PLN/pln_eval/out/full_comparison_table_v2.csv

bmr_dir  <- "/projects/genomic-ml/da2343/PLN/pln_eval/bmr"
data_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval/data"
out_dir  <- "/projects/genomic-ml/da2343/PLN/pln_eval/out"

# ── Helpers ──────────────────────────────────────────────────────────────────

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

compute_mac <- function(counts) {
  N      <- nrow(counts)
  p      <- ncol(counts)
  X      <- log1p(counts)
  cm     <- cor(X)
  diag(cm) <- 0
  mac    <- mean(abs(cm), na.rm = TRUE)
  data.frame(N = N, p = p, N_over_p = round(N / p, 2),
             mac = round(mac, 3), stringsAsFactors = FALSE)
}

# ── Mapping: bmr CSV name → data file path (for new datasets only) ───────────
# Datasets already in covariance_metrics.csv are handled by the name map below.

new_dataset_paths <- list(
  "amgut1" =
    file.path(data_dir, "amgut1_update.csv"),
  "amgut2" =
    file.path(data_dir, "amgut2_update.csv"),
  "hiv_lozupone_family" =
    file.path(data_dir, "family/microbiomehd/hiv_lozupone_family.tsv.gz"),
  "mbqc_integrated_otus_family" =
    file.path(data_dir, "family/mbqc/mbqc_integrated_otus_family.tsv.gz"),
  "crc_zhao_genus" =
    file.path(data_dir, "genus/microbiomehd/crc_zhao_genus.tsv.gz"),
  "ibd_gevers_2014_family" =
    file.path(data_dir, "family/microbiomehd/ibd_gevers_2014_family.tsv.gz"),
  "nash_chan_family" =
    file.path(data_dir, "family/microbiomehd/nash_chan_family.tsv.gz"),
  "ob_ross_genus" =
    file.path(data_dir, "genus/microbiomehd/ob_ross_genus.tsv.gz"),
  "par_scheperjans_family" =
    file.path(data_dir, "family/microbiomehd/par_scheperjans_family.tsv.gz"),
  "t1d_alkanani_family" =
    file.path(data_dir, "family/microbiomehd/t1d_alkanani_family.tsv.gz"),
  "cdi_schubert_family" =
    file.path(data_dir, "family/microbiomehd/cdi_schubert_family.tsv.gz"),
  "20170526.Castro-NallarE_2015.metaphlan_bugs_list.oralcavity_genus" =
    file.path(data_dir, "genus/curated_metagenomic_data/20170526.Castro-NallarE_2015.metaphlan_bugs_list.oralcavity_genus.tsv.gz"),
  "20170526.ChngKR_2016.metaphlan_bugs_list.skin_genus" =
    file.path(data_dir, "genus/curated_metagenomic_data/20170526.ChngKR_2016.metaphlan_bugs_list.skin_genus.tsv.gz"),
  "20181025.OlmMR_2017.metaphlan_bugs_list.stool_genus" =
    file.path(data_dir, "genus/curated_metagenomic_data/20181025.OlmMR_2017.metaphlan_bugs_list.stool_genus.tsv.gz"),
  "2021-10-14.MehtaRS_2018.relative_abundance_genus" =
    file.path(data_dir, "genus/curated_metagenomic_data/2021-10-14.MehtaRS_2018.relative_abundance_genus.tsv.gz")
)

# ── Mapping: bmr CSV name → covariance_metrics dataset name ──────────────────

cov_name_map <- c(
  "2021-10-14.AsnicarF_2017.relative_abundance_species"        = "AsnicarF_2017 species",
  "20170526.HMP_2012.metaphlan_bugs_list.nasalcavity_species"  = "HMP_2012 nasalcavity species",
  "HMP_2012.metaphlan_bugs_list.anterior_nares_species"        = "HMP_2012 anterior_nares species",
  "2022-04-13.PetersBA_2019.relative_abundance_genus"          = "PetersBA_2019 genus",
  "2021-03-31.OlmMR_2017.relative_abundance_genus"             = "OlmMR_2017 genus",
  "20180425.KarlssonFH_2013.metaphlan_bugs_list.stool_genus"   = "KarlssonFH_2013 stool genus",
  "20181025.YuJ_2015.metaphlan_bugs_list.stool_genus"          = "YuJ_2015 stool genus",
  "20190424.CosteaPI_2017.metaphlan_bugs_list.stool_genus"     = "CosteaPI_2017 stool genus",
  "2021-03-31.MehtaRS_2018.relative_abundance_genus"           = "MehtaRS_2018 genus",
  "2021-03-31.ShaoY_2019.relative_abundance_genus"             = "ShaoY_2019 genus",
  "2021-03-31.VilaAV_2018.relative_abundance_genus"            = "VilaAV_2018 genus",
  "2021-10-14.RubelMA_2020.relative_abundance_species"         = "RubelMA_2020 species",
  "ra_littman_family"                                           = "ra_littman family",
  "otu_table_psn_v13_family"                                    = "otu_table_psn_v13 family",
  "ob_zeevi_family"                                             = "ob_zeevi family",
  "nash_ob_baker_family"                                        = "nash_ob_baker family",
  "diabimmune_karelia_16s_family"                               = "diabimmune_karelia family",
  "crc_baxter"                                                  = "crc_baxter",
  "crc_zeller"                                                  = "crc_zeller",
  "cdi_schubert_genus"                                          = "cdi_schubert genus",
  "BertuzziAS_2018_genus"                                       = "BertuzziAS_2018 genus"
)

# ── Load covariance metrics ───────────────────────────────────────────────────

cov <- read.csv(file.path(out_dir, "covariance_metrics.csv"),
                stringsAsFactors = FALSE)

# ── Find all bmr CSVs (exclude sub-sampled and collect_amgut2) ───────────────

all_csvs <- list.files(bmr_dir, pattern = "\\.csv$", full.names = FALSE)
all_csvs <- all_csvs[!grepl("_[0-9]+_samples\\.csv$|collect_amgut2", all_csvs)]
dataset_names <- sub("\\.csv$", "", all_csvs)

# ── Process each dataset ──────────────────────────────────────────────────────

results <- list()

for (ds in dataset_names) {
  cat("Processing:", ds, "...")

  # Read benchmark results
  bmr_path <- file.path(bmr_dir, paste0(ds, ".csv"))
  bmr      <- read.csv(bmr_path, stringsAsFactors = FALSE)

  featureless <- bmr$mean_deviance[bmr$learner_id == "regr.featureless"]
  lasso       <- bmr$mean_deviance[bmr$learner_id == "regr.lasso"]
  pln         <- bmr$mean_deviance[bmr$learner_id == "regr.pln"]

  if (length(featureless) == 0 || length(lasso) == 0 || length(pln) == 0) {
    cat(" SKIP (missing learners)\n")
    next
  }

  # Get or compute N, p, MAC
  if (ds %in% names(cov_name_map)) {
    cov_name <- cov_name_map[ds]
    row      <- cov[cov$dataset == cov_name, ]
    if (nrow(row) == 0) {
      cat(" SKIP (not found in covariance_metrics)\n")
      next
    }
    N       <- row$N
    p       <- row$p
    N_over_p <- row$N_over_p
    mac     <- row$mac
  } else if (ds %in% names(new_dataset_paths)) {
    dpath <- new_dataset_paths[[ds]]
    if (!file.exists(dpath)) {
      cat(" SKIP (data file not found:", dpath, ")\n")
      next
    }
    if (grepl("\\.csv$", dpath)) {
      counts <- load_counts_csv(dpath)
    } else {
      counts <- load_counts_tsv(dpath)
    }
    m       <- compute_mac(counts)
    N       <- m$N
    p       <- m$p
    N_over_p <- m$N_over_p
    mac     <- m$mac
  } else {
    cat(" SKIP (no data source mapped)\n")
    next
  }

  winner           <- if (pln < lasso) "PLN" else "Lasso"
  pln_vs_lasso_pct <- round((lasso - pln) / lasso * 100, 1)

  results[[length(results) + 1]] <- data.frame(
    dataset          = ds,
    N                = N,
    p                = p,
    N_over_p         = N_over_p,
    winner           = winner,
    pln_vs_lasso_pct = paste0(ifelse(pln_vs_lasso_pct > 0, "+", ""), pln_vs_lasso_pct, "%"),
    featureless      = round(featureless, 4),
    lasso            = round(lasso, 4),
    pln              = round(pln, 4),
    mac              = mac,
    stringsAsFactors = FALSE
  )
  cat(" done (N=", N, "p=", p, "mac=", mac, "winner=", winner, ")\n")
}

out <- do.call(rbind, results)
out <- out[order(out$N_over_p), ]

out_path <- file.path(out_dir, "full_comparison_table_v2.csv")
write.csv(out, out_path, row.names = FALSE)
cat("\nTotal datasets:", nrow(out), "\n")
cat("PLN wins:", sum(out$winner == "PLN"), "\n")
cat("Lasso wins:", sum(out$winner == "Lasso"), "\n")
cat("Saved to:", out_path, "\n")
