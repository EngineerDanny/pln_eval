#!/usr/bin/env Rscript
# bootstrap_f1.R
#
# For each of 5 network-truth datasets, runs B=20 bootstrap resamples of the
# abundance matrix rows, fits PLNnetwork and LOTO_glmnet_CV1se on each
# resample, and computes F1 against the FIXED ground truth.
#
# Output: /projects/genomic-ml/da2343/PLN/pln_eval/out/bootstrap_f1_results.csv
#
# Run via SLURM:
#   srun /projects/genomic-ml/da2343/soak-r/bin/Rscript bootstrap_f1.R

suppressPackageStartupMessages({
  library(data.table)
  library(PLNmodels)
  library(glassoFast)
  library(glmnet)
})

set.seed(42)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
base_dir  <- "/projects/genomic-ml/da2343/PLN/pln_eval"
truth_dir <- file.path(base_dir, "data/interaction_ground_truth")
out_dir   <- file.path(base_dir, "out")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Datasets and bootstrap settings
# ---------------------------------------------------------------------------
datasets <- c(
  "omm12",
  "omm12_keystone_2023",
  "pairinterax",
  "butyrate_assembly_2021",
  "host_fitness_2018"
)

B <- 20L   # number of bootstrap resamples per dataset × method

# ---------------------------------------------------------------------------
# Helper utilities (copied verbatim from fit_truth_benchmark.R)
# ---------------------------------------------------------------------------

read_tsv_gz <- function(path) {
  read.delim(gzfile(path), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

sort_pair <- function(a, b) {
  ifelse(a <= b, paste(a, b, sep = "||"), paste(b, a, sep = "||"))
}

prepare_count_like <- function(X, target_depth = 10000L) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  X[!is.finite(X)] <- 0
  X[X < 0] <- 0

  is_integerish <- all(abs(X - round(X)) < 1e-8)
  if (is_integerish) {
    return(round(X))
  }

  rs <- rowSums(X)
  scaled <- matrix(0, nrow(X), ncol(X), dimnames = dimnames(X))
  keep <- rs > 0
  scaled[keep, ] <- round(X[keep, , drop = FALSE] / rs[keep] * target_depth)
  scaled
}

pln_prepare_counts <- function(counts) {
  prepare_data(
    counts = counts,
    covariates = data.frame(Intercept = rep(1, nrow(counts)), row.names = rownames(counts)),
    offset = "none"
  )
}

loto_to_omega <- function(coef_matrix, sigma_diag) {
  p <- nrow(coef_matrix)
  Omega <- matrix(0, p, p)
  for (j in seq_len(p)) {
    Omega[j, j] <- 1 / sigma_diag[j]
    for (k in which(coef_matrix[j, ] != 0)) {
      Omega[j, k] <- -coef_matrix[j, k] * Omega[j, j]
    }
  }
  Omega <- (Omega + t(Omega)) / 2
  eig <- eigen(Omega, symmetric = TRUE)
  eig$values <- pmax(eig$values, 1e-6)
  eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
}

support_from_coef_mat <- function(coef_mat, require_mutual = TRUE, tol = 1e-8) {
  nz <- abs(coef_mat) > tol
  diag(nz) <- FALSE
  support <- if (require_mutual) nz & t(nz) else nz | t(nz)
  diag(support) <- FALSE
  support
}

omega_to_prediction <- function(Omega, taxa) {
  idx <- which(upper.tri(Omega) & Omega != 0, arr.ind = TRUE)
  if (nrow(idx) == 0L) {
    return(data.frame(taxon_1 = character(), taxon_2 = character(),
                      sign = integer(), weight = numeric(),
                      stringsAsFactors = FALSE))
  }
  w <- Omega[idx]
  data.frame(
    taxon_1 = taxa[idx[, 1]],
    taxon_2 = taxa[idx[, 2]],
    sign    = ifelse(w > 0, 1L, -1L),
    weight  = as.numeric(w),
    stringsAsFactors = FALSE
  )
}

support_to_prediction <- function(support, weight_mat, taxa) {
  idx <- which(upper.tri(support) & support, arr.ind = TRUE)
  if (nrow(idx) == 0L) {
    return(data.frame(taxon_1 = character(), taxon_2 = character(),
                      sign = integer(), weight = numeric(),
                      stringsAsFactors = FALSE))
  }
  w <- weight_mat[idx]
  data.frame(
    taxon_1 = taxa[idx[, 1]],
    taxon_2 = taxa[idx[, 2]],
    sign    = ifelse(w > 0, 1L, -1L),
    weight  = as.numeric(w),
    stringsAsFactors = FALSE
  )
}

compute_metrics <- function(truth, pred, dataset_name) {
  truth$key <- sort_pair(truth$taxon_1, truth$taxon_2)
  pred$key  <- sort_pair(pred$taxon_1,  pred$taxon_2)
  truth_keys <- unique(truth$key)
  pred_keys  <- unique(pred$key)
  tp_keys    <- intersect(truth_keys, pred_keys)
  fp_keys    <- setdiff(pred_keys, truth_keys)
  fn_keys    <- setdiff(truth_keys, pred_keys)
  tp <- length(tp_keys); fp <- length(fp_keys); fn <- length(fn_keys)
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  recall    <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  f1        <- if (is.finite(precision + recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }
  data.frame(
    dataset   = dataset_name,
    n_truth   = length(truth_keys),
    n_pred    = length(pred_keys),
    tp = tp, fp = fp, fn = fn,
    precision = precision,
    recall    = recall,
    f1        = f1,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Fitting functions (copied verbatim from fit_truth_benchmark.R)
# ---------------------------------------------------------------------------

fit_plnnetwork <- function(X_train_raw) {
  pln_init <- PLN(
    Abundance ~ 1,
    data = pln_prepare_counts(X_train_raw),
    control = PLN_param(backend = "nlopt", trace = 0)
  )
  sigma_pln <- pln_init$model_par$Sigma
  max_rho   <- max(abs(sigma_pln[upper.tri(sigma_pln)]))
  rho_seq   <- exp(seq(log(max_rho), log(max_rho * 0.005), length.out = 30))

  fit <- PLNnetwork(
    Abundance ~ 1,
    data      = pln_prepare_counts(X_train_raw),
    penalties = rho_seq,
    control   = PLNnetwork_param(
      backend           = "nlopt",
      trace             = 0,
      penalize_diagonal = FALSE,
      inception         = pln_init
    )
  )
  fit$getBestModel(crit = "EBIC")$model_par$Omega
}

safe_glmnet_coefs <- function(x_train, y_train) {
  p_minus_1 <- ncol(x_train)
  zero_coefs <- rep(0, p_minus_1)
  if (p_minus_1 == 0L || length(unique(y_train)) < 2L) return(zero_coefs)
  x_train   <- as.matrix(x_train)
  cv_nfolds <- max(3L, min(5L, nrow(x_train)))
  fit <- tryCatch(
    cv.glmnet(x_train, y_train, family = "poisson", alpha = 1, nfolds = cv_nfolds),
    error = function(e) NULL
  )
  if (is.null(fit)) return(zero_coefs)
  coefs <- tryCatch(as.numeric(coef(fit, s = "lambda.1se")[-1]), error = function(e) zero_coefs)
  if (length(coefs) != p_minus_1) zero_coefs else coefs
}

fit_loto_glmnet <- function(X_train_raw, X_train_log) {
  p        <- ncol(X_train_raw)
  coef_mat <- matrix(0, p, p)
  for (j in seq_len(p)) {
    y_train <- X_train_raw[, j]
    x_train <- X_train_log[, -j, drop = FALSE]
    coef_mat[j, -j] <- safe_glmnet_coefs(x_train, y_train)
  }
  sigma_diag <- pmax(apply(X_train_log, 2, var), 1e-6)
  list(omega = loto_to_omega(coef_mat, sigma_diag), coef_mat = coef_mat)
}

# ---------------------------------------------------------------------------
# Data loading (mirrors load_benchmark_matrix() from fit_truth_benchmark.R,
# but scoped to work with a dataset name rather than a global proc_dir)
# ---------------------------------------------------------------------------

load_dataset <- function(dataset) {
  proc_dir <- file.path(truth_dir, dataset, "processed")

  # Abundance matrix
  if (dataset %in% c("omm12", "omm12_keystone_2023")) {
    df <- read_tsv_gz(file.path(proc_dir, "community_absabundance_in_vitro.tsv.gz"))
  } else {
    df <- read_tsv_gz(file.path(proc_dir, "abundance_matrix.tsv.gz"))
  }
  X <- as.matrix(df[, -1, drop = FALSE])
  rownames(X) <- df$sample_id

  # Fixed ground truth (never resampled)
  truth <- read_tsv_gz(file.path(proc_dir, "truth_undirected.tsv.gz"))

  list(X = X, truth = truth)
}

# ---------------------------------------------------------------------------
# Main bootstrap loop
# ---------------------------------------------------------------------------

all_rows <- list()

for (dataset in datasets) {
  message("\n=== Dataset: ", dataset, " ===")

  dat   <- load_dataset(dataset)
  X_raw <- dat$X
  truth <- dat$truth

  # Pre-process the FULL matrix once to get the column set and filter
  X_counts_full <- prepare_count_like(X_raw)
  col_keep      <- colSums(X_counts_full) > 0
  X_counts_full <- X_counts_full[rowSums(X_counts_full) > 0, col_keep, drop = FALSE]
  taxa          <- colnames(X_counts_full)
  n_full        <- nrow(X_counts_full)

  for (b in seq_len(B)) {
    message("  Bootstrap rep ", b, " / ", B)

    # Resample rows with replacement
    boot_idx <- sample(n_full, replace = TRUE)
    X_boot   <- X_counts_full[boot_idx, , drop = FALSE]

    # Drop any rows that became all-zero after resampling (extremely rare)
    row_keep <- rowSums(X_boot) > 0
    X_boot   <- X_boot[row_keep, , drop = FALSE]

    # Unique row names to satisfy PLNmodels
    rownames(X_boot) <- paste0("bs_", seq_len(nrow(X_boot)))

    X_log_boot <- log1p(X_boot)

    # --- PLNnetwork ---
    f1_pln <- tryCatch({
      omega <- fit_plnnetwork(X_boot)
      pred  <- omega_to_prediction(omega, taxa)
      m     <- compute_metrics(truth, pred, dataset)
      m$f1
    }, error = function(e) {
      message("    PLNnetwork error: ", conditionMessage(e))
      NA_real_
    })

    # --- LOTO_glmnet_CV1se ---
    f1_glmnet <- tryCatch({
      fit     <- fit_loto_glmnet(X_boot, X_log_boot)
      support <- support_from_coef_mat(fit$coef_mat, require_mutual = TRUE)
      wmat    <- (fit$coef_mat + t(fit$coef_mat)) / 2
      pred    <- support_to_prediction(support, wmat, taxa)
      m       <- compute_metrics(truth, pred, dataset)
      m$f1
    }, error = function(e) {
      message("    LOTO_glmnet_CV1se error: ", conditionMessage(e))
      NA_real_
    })

    all_rows <- c(all_rows, list(
      data.frame(dataset = dataset, method = "PLNnetwork",        bootstrap_rep = b, f1 = f1_pln,    stringsAsFactors = FALSE),
      data.frame(dataset = dataset, method = "LOTO_glmnet_CV1se", bootstrap_rep = b, f1 = f1_glmnet, stringsAsFactors = FALSE)
    ))
  }
}

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------

results <- rbindlist(all_rows, fill = FALSE)
out_path <- file.path(out_dir, "bootstrap_f1_results.csv")
fwrite(results, out_path)
message("\nSaved results to: ", out_path)

# ---------------------------------------------------------------------------
# Summary table: dataset × method → mean_f1, se_f1
# ---------------------------------------------------------------------------

summary_dt <- results[, .(
  mean_f1 = mean(f1, na.rm = TRUE),
  se_f1   = sd(f1,   na.rm = TRUE) / sqrt(sum(!is.na(f1)))
), keyby = .(dataset, method)]

message("\n=== Bootstrap F1 Summary (B=", B, ") ===")
print(summary_dt)
