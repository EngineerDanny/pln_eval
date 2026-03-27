library(data.table)
library(PLNmodels)
library(glassoFast)
library(glmnet)

args <- commandArgs(trailingOnly = TRUE)
base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
data_dir <- file.path(base_dir, "data")

load_counts <- function(dataname_or_path, data_dir) {
  if (grepl("\\.tsv\\.gz$", dataname_or_path)) {
    if (file.exists(dataname_or_path)) {
      tsv_path <- dataname_or_path
    } else {
      tsv_path <- file.path(data_dir, dataname_or_path)
    }
    if (!file.exists(tsv_path)) stop("File not found: ", dataname_or_path)
    raw <- read.table(
      gzfile(tsv_path),
      sep = "\t",
      header = TRUE,
      row.names = 1,
      check.names = FALSE
    )
    counts <- t(as.matrix(raw))
    return(counts[rowSums(counts) > 0, colSums(counts) > 0, drop = FALSE])
  }

  if (grepl("\\.csv$", dataname_or_path)) {
    if (file.exists(dataname_or_path)) {
      csv_path <- dataname_or_path
    } else {
      csv_path <- file.path(data_dir, dataname_or_path)
    }
    if (!file.exists(csv_path)) stop("File not found: ", dataname_or_path)
    dt <- fread(csv_path)
    if ("Group_ID" %in% names(dt)) dt[, Group_ID := NULL]
    return(as.matrix(dt))
  }

  tsv_path <- list.files(
    data_dir,
    pattern = paste0(dataname_or_path, ".*\\.tsv\\.gz$"),
    recursive = TRUE,
    full.names = TRUE
  )
  csv_path <- file.path(data_dir, paste0(dataname_or_path, "_update.csv"))

  if (length(tsv_path) > 0) {
    raw <- read.table(
      gzfile(tsv_path[1]),
      sep = "\t",
      header = TRUE,
      row.names = 1,
      check.names = FALSE
    )
    counts <- t(as.matrix(raw))
    return(counts[rowSums(counts) > 0, colSums(counts) > 0, drop = FALSE])
  }
  if (file.exists(csv_path)) {
    dt <- fread(csv_path)
    if ("Group_ID" %in% names(dt)) dt[, Group_ID := NULL]
    return(as.matrix(dt))
  }

  stop("No data file found for: ", dataname_or_path)
}

dataset_arg <- if (length(args) >= 1) args[1] else "amgut2_update.csv"
dataset_tag <- if (grepl("\\.tsv\\.gz$", dataset_arg)) {
  gsub("\\.tsv\\.gz$", "", basename(dataset_arg))
} else if (grepl("\\.csv$", dataset_arg)) {
  gsub("\\.csv$", "", basename(dataset_arg))
} else {
  dataset_arg
}

out_dir <- file.path(base_dir, "out", "network_graph", dataset_tag)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pln_prepare_counts <- function(counts) {
  prepare_data(
    counts = counts,
    covariates = data.frame(
      Intercept = rep(1, nrow(counts)),
      row.names = rownames(counts)
    ),
    offset = "none"
  )
}

omega_to_edge_dt <- function(Omega, taxa, rep_id, method) {
  idx <- which(upper.tri(Omega) & Omega != 0, arr.ind = TRUE)
  if (nrow(idx) == 0L) {
    return(data.table(
      rep = integer(),
      method = character(),
      taxon_i = character(),
      taxon_j = character(),
      i = integer(),
      j = integer(),
      weight = numeric(),
      abs_weight = numeric(),
      sign = character()
    ))
  }

  weights <- Omega[idx]
  data.table(
    rep = rep_id,
    method = method,
    taxon_i = taxa[idx[, 1]],
    taxon_j = taxa[idx[, 2]],
    i = idx[, 1],
    j = idx[, 2],
    weight = weights,
    abs_weight = abs(weights),
    sign = ifelse(weights > 0, "positive", "negative")
  )[order(-abs_weight)]
}

fit_pln_glasso <- function(X_train_raw) {
  n <- nrow(X_train_raw)
  final_model <- PLN(
    Abundance ~ 1,
    data = pln_prepare_counts(X_train_raw),
    control = PLN_param(backend = "nlopt", trace = 0)
  )
  sigma_pln <- final_model$model_par$Sigma
  p <- ncol(sigma_pln)

  max_rho <- max(abs(sigma_pln[upper.tri(sigma_pln)]))
  rho_seq <- exp(seq(log(max_rho), log(max_rho * 0.005), length.out = 30))
  best_rho <- rho_seq[1]
  best_ebic <- Inf

  for (rho in rho_seq) {
    tryCatch({
      gl <- glassoFast(sigma_pln, rho = rho)
      logdet <- determinant(gl$wi, logarithm = TRUE)$modulus
      loglik_train <- (n / 2) * (logdet - sum(diag(gl$wi %*% sigma_pln)))
      n_edges <- sum(gl$wi[upper.tri(gl$wi)] != 0)
      ebic <- -2 * loglik_train + n_edges * log(n) + 4 * 0.5 * n_edges * log(p)
      if (ebic < best_ebic) {
        best_ebic <- ebic
        best_rho <- rho
      }
    }, error = function(e) NULL)
  }

  gl_best <- glassoFast(sigma_pln, rho = best_rho)
  list(
    omega = gl_best$wi,
    best_rho = best_rho
  )
}

fit_loto_pln_glasso <- function(X_train_raw, X_train_log, rho_grid_len = 12L) {
  n <- nrow(X_train_raw)
  p <- ncol(X_train_raw)
  coef_mat <- matrix(0, p, p)
  sigma_diag <- pmax(apply(X_train_log, 2, var), 1e-6)

  for (target_idx in seq_len(p)) {
    order_idx <- c(target_idx, setdiff(seq_len(p), target_idx))
    counts_target <- X_train_raw[, order_idx, drop = FALSE]
    log_target <- X_train_log[, order_idx, drop = FALSE]

    fit <- PLN(
      Abundance ~ 1,
      data = pln_prepare_counts(counts_target),
      control = PLN_param(backend = "nlopt", trace = 0)
    )
    sigma_pln <- fit$model_par$Sigma
    max_rho <- max(abs(sigma_pln[upper.tri(sigma_pln)]))
    if (!is.finite(max_rho) || max_rho <= 0) next

    rho_seq <- exp(seq(log(max_rho), log(max_rho * 0.01), length.out = rho_grid_len))
    centered <- sweep(log_target, 2, colMeans(log_target))
    best_bic <- Inf
    best_omega <- NULL

    for (rho in rho_seq) {
      tryCatch({
        gl <- glassoFast(sigma_pln, rho = rho)
        omega <- gl$wi
        omega_jj <- omega[1, 1]
        if (!is.finite(omega_jj) || omega_jj <= 0) return(NULL)

        cond_mean <- -centered[, -1, drop = FALSE] %*% omega[-1, 1] / omega_jj
        resid <- centered[, 1] - as.vector(cond_mean)
        cond_ll_total <- sum(
          0.5 * log(omega_jj) - 0.5 * omega_jj * resid^2 - 0.5 * log(2 * pi)
        )
        df <- sum(omega[1, -1] != 0)
        bic <- -2 * cond_ll_total + df * log(n)

        if (bic < best_bic) {
          best_bic <- bic
          best_omega <- omega
        }
      }, error = function(e) NULL)
    }

    if (is.null(best_omega)) next

    coef_mat[target_idx, setdiff(seq_len(p), target_idx)] <- -best_omega[1, -1] / best_omega[1, 1]
    sigma_diag[target_idx] <- 1 / best_omega[1, 1]
  }

  list(
    omega = loto_to_omega(coef_mat, sigma_diag),
    coef_mat = coef_mat
  )
}

## ── Setup ──────────────────────────────────────────────────────
counts_full <- load_counts(dataset_arg, data_dir)
col_sums <- colSums(counts_full)
top_taxa <- order(col_sums, decreasing = TRUE)
X_raw <- counts_full[, top_taxa]
X_log <- log1p(X_raw)
rownames(X_raw) <- paste0("S", seq_len(nrow(X_raw)))
rownames(X_log) <- rownames(X_raw)

n_total <- nrow(X_log)
p <- ncol(X_log)
cat(sprintf("Data: %s | %d samples, %d taxa (log1p)\n", dataset_tag, n_total, p))

## ── Parameters ─────────────────────────────────────────────────
n_reps     <- 2
n_train    <- 200
set.seed(42)

## ── Held-out pseudo-likelihood under multivariate Gaussian ─────
held_out_pseudologlik <- function(Omega, mu, X_test) {
  X_centered <- sweep(X_test, 2, mu)
  p <- ncol(X_test)
  cond_ll <- numeric(p)

  for (j in seq_len(p)) {
    omega_jj <- Omega[j, j]
    if (!is.finite(omega_jj) || omega_jj <= 0) return(NA_real_)

    cond_mean <- -X_centered[, -j, drop = FALSE] %*% Omega[-j, j] / omega_jj
    resid <- X_centered[, j] - as.vector(cond_mean)
    cond_ll[j] <- mean(
      0.5 * log(omega_jj) - 0.5 * omega_jj * resid^2 - 0.5 * log(2 * pi)
    )
  }

  mean(cond_ll)
}

## ── LOTO: reconstruct precision matrix from neighborhood coefs ─
# For each taxon j, coef vector beta_j gives regression on others.
# Reconstruct Omega using: Omega[j,-j] = -beta_j * Omega[j,j]
# We use a simple symmetrized approach from the AND-rule adjacency.
loto_to_omega <- function(coef_matrix, sigma_diag) {
  p <- nrow(coef_matrix)
  Omega <- matrix(0, p, p)
  for (j in 1:p) {
    Omega[j, j] <- 1 / sigma_diag[j]
    for (k in which(coef_matrix[j, ] != 0)) {
      Omega[j, k] <- -coef_matrix[j, k] * Omega[j, j]
    }
  }
  # Symmetrize
  Omega <- (Omega + t(Omega)) / 2
  # Ensure positive definite
  eig <- eigen(Omega, symmetric = TRUE)
  eig$values <- pmax(eig$values, 1e-6)
  Omega_pd <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  Omega_pd
}

## ── Run replicates ─────────────────────────────────────────────
results_list <- list()
edge_results <- list()

for (r in 1:n_reps) {
  cat(sprintf("\n=== Replicate %d ===\n", r))
  train_idx <- sample(n_total, n_train)
  test_idx  <- setdiff(1:n_total, train_idx)

  X_train_raw <- X_raw[train_idx, ]
  X_train <- X_log[train_idx, ]
  X_test  <- X_log[test_idx, ]
  mu_train <- colMeans(X_train)

  ## ── Method 1: PLN + glasso ───────────────────────────────────
  cat("  PLN+glasso... ")
  t1 <- proc.time()
  fit1 <- fit_pln_glasso(X_train_raw)
  omega1 <- fit1$omega
  ll1 <- held_out_pseudologlik(omega1, mu_train, X_test)
  n_edges1 <- sum(omega1[upper.tri(omega1)] != 0)
  t1e <- (proc.time() - t1)["elapsed"]
  cat(sprintf(
    "loglik=%.2f, edges=%d, rho=%.4f, time=%.1fs\n",
    ll1, n_edges1, fit1$best_rho, t1e
  ))

  ## ── Method 2: LOTO PLN + glasso ──────────────────────────────
  cat("  LOTO PLN+glasso... ")
  t2 <- proc.time()
  fit2 <- fit_loto_pln_glasso(X_train_raw, X_train)
  omega2 <- fit2$omega
  ll2 <- held_out_pseudologlik(omega2, mu_train, X_test)
  n_edges2 <- sum((fit2$coef_mat != 0 & t(fit2$coef_mat) != 0)[upper.tri(fit2$coef_mat)])
  t2e <- (proc.time() - t2)["elapsed"]
  cat(sprintf("loglik=%.2f, edges=%d, time=%.1fs\n", ll2, n_edges2, t2e))

  ## ── Method 3: LOTO glmnet + CV1se ────────────────────────────
  cat("  LOTO glmnet+CV1se... ")
  t3 <- proc.time()
  coef_mat <- matrix(0, p, p)
  for (j in 1:p) {
    y_train <- X_raw[train_idx, j]
    x_train <- X_train[, -j]
    cvfit <- cv.glmnet(x_train, y_train, family = "poisson", alpha = 1)
    coefs <- as.numeric(coef(cvfit, s = "lambda.1se")[-1])
    coef_mat[j, -j] <- coefs
  }
  sigma_diag <- apply(X_train, 2, var)
  omega3 <- loto_to_omega(coef_mat, sigma_diag)
  ll3 <- held_out_pseudologlik(omega3, mu_train, X_test)
  n_edges3 <- sum((coef_mat != 0 & t(coef_mat) != 0)[upper.tri(coef_mat)])
  t3e <- (proc.time() - t3)["elapsed"]
  cat(sprintf("loglik=%.2f, edges=%d, time=%.1fs\n", ll3, n_edges3, t3e))

  ## ── Baseline: diagonal Omega (independence, no edges) ─────────
  omega0_diag <- diag(1 / apply(X_train, 2, var))
  ll0_diag <- held_out_pseudologlik(omega0_diag, mu_train, X_test)

  omega_map <- list(
    Baseline_diagonal = omega0_diag,
    PLN_glasso = omega1,
    LOTO_PLN_glasso = omega2,
    LOTO_glmnet_CV1se = omega3
  )

  for (method_name in names(omega_map)) {
    edge_results[[length(edge_results) + 1L]] <- omega_to_edge_dt(
      omega_map[[method_name]],
      taxa = colnames(X_train),
      rep_id = r,
      method = method_name
    )
  }

  results_list[[r]] <- data.table(
    rep    = r,
    method = c(
      "Baseline_diagonal",
      "PLN_glasso",
      "LOTO_PLN_glasso",
      "LOTO_glmnet_CV1se"
    ),
    loglik  = round(c(ll0_diag, ll1, ll2, ll3), 3),
    n_edges = c(0, n_edges1, n_edges2, n_edges3)
  )
}

## ── Summary ────────────────────────────────────────────────────
all_results <- rbindlist(results_list)
edge_dt <- rbindlist(edge_results, fill = TRUE)
summary_dt <- all_results[, .(
  mean_loglik = round(mean(loglik), 3),
  sd_loglik   = round(sd(loglik), 3),
  mean_edges  = round(mean(n_edges), 1)
), by = method][order(-mean_loglik)]

fwrite(
  all_results,
  file.path(out_dir, "network_graph_pseudologlik_replicates.csv")
)
fwrite(
  summary_dt,
  file.path(out_dir, "network_graph_pseudologlik_summary.csv")
)
fwrite(
  edge_dt,
  file.path(out_dir, "network_graph_edges.csv")
)

cat("\n=== Per-replicate results ===\n")
print(dcast(all_results, method ~ rep, value.var = "loglik"))

cat("\n=== Summary (higher pseudo-loglik = better) ===\n")
print(summary_dt[order(-mean_loglik)])
cat("\nSaved:\n")
cat(file.path(out_dir, "network_graph_pseudologlik_replicates.csv"), "\n")
cat(file.path(out_dir, "network_graph_pseudologlik_summary.csv"), "\n")
cat(file.path(out_dir, "network_graph_edges.csv"), "\n")
cat("\nDone.\n")
