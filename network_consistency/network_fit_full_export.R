library(data.table)
library(PLNmodels)
library(glassoFast)
library(glmnet)

args <- commandArgs(trailingOnly = TRUE)
dataset_arg <- if (length(args) >= 1) args[1] else "amgut2_update.csv"
top_n_taxa <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_

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

dataset_tag <- if (grepl("\\.tsv\\.gz$", dataset_arg)) {
  gsub("\\.tsv\\.gz$", "", basename(dataset_arg))
} else if (grepl("\\.csv$", dataset_arg)) {
  gsub("\\.csv$", "", basename(dataset_arg))
} else {
  dataset_arg
}

out_tag <- if (length(args) >= 3) {
  args[3]
} else if (is.na(top_n_taxa)) {
  sprintf("%s_full_fit", dataset_tag)
} else {
  sprintf("%s_top_%d", dataset_tag, top_n_taxa)
}

out_dir <- file.path(base_dir, "out", "network_graph", out_tag)
omega_dir <- file.path(out_dir, "omegas")
dir.create(omega_dir, recursive = TRUE, showWarnings = FALSE)

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

loto_to_omega <- function(coef_matrix, sigma_diag) {
  p <- nrow(coef_matrix)
  Omega <- matrix(0, p, p)
  for (j in 1:p) {
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

safe_glmnet_coefs <- function(x_train, y_train) {
  p_minus_1 <- ncol(x_train)
  zero_coefs <- rep(0, p_minus_1)
  if (p_minus_1 == 0L) return(zero_coefs)
  if (length(unique(y_train)) < 2L) return(zero_coefs)

  x_train <- as.matrix(x_train)
  cv_nfolds <- max(3L, min(5L, nrow(x_train)))

  fit <- tryCatch(
    cv.glmnet(
      x_train,
      y_train,
      family = "poisson",
      alpha = 1,
      nfolds = cv_nfolds
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(zero_coefs)

  coefs <- tryCatch(
    as.numeric(coef(fit, s = "lambda.1se")[-1]),
    error = function(e) zero_coefs
  )

  if (length(coefs) != p_minus_1) return(zero_coefs)
  coefs
}

fit_loto_glmnet <- function(X_train_raw, X_train_log) {
  p <- ncol(X_train_raw)
  coef_mat <- matrix(0, p, p)
  for (j in seq_len(p)) {
    y_train <- X_train_raw[, j]
    x_train <- X_train_log[, -j, drop = FALSE]
    coef_mat[j, -j] <- safe_glmnet_coefs(x_train, y_train)
  }
  sigma_diag <- apply(X_train_log, 2, var)
  list(
    omega = loto_to_omega(coef_mat, sigma_diag),
    coef_mat = coef_mat
  )
}

omega_to_edge_dt <- function(Omega, taxa, method) {
  idx <- which(upper.tri(Omega) & Omega != 0, arr.ind = TRUE)
  if (nrow(idx) == 0L) {
    return(data.table(
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

counts_full <- load_counts(dataset_arg, data_dir)
col_sums <- colSums(counts_full)
taxa_order <- order(col_sums, decreasing = TRUE)
if (!is.na(top_n_taxa)) {
  taxa_order <- taxa_order[seq_len(min(top_n_taxa, length(taxa_order)))]
}
X_raw <- counts_full[, taxa_order, drop = FALSE]
X_log <- log1p(X_raw)
rownames(X_raw) <- paste0("S", seq_len(nrow(X_raw)))
rownames(X_log) <- rownames(X_raw)

cat(sprintf(
  "Full-data fit: %s | %d samples, %d taxa (tag: %s)\n",
  dataset_tag, nrow(X_raw), ncol(X_raw), out_tag
))

fit_pln <- fit_pln_glasso(X_raw)
fit_loto_pln <- fit_loto_pln_glasso(X_raw, X_log)
fit_loto_glm <- fit_loto_glmnet(X_raw, X_log)

omega_map <- list(
  PLN_glasso = fit_pln$omega,
  LOTO_PLN_glasso = fit_loto_pln$omega,
  LOTO_glmnet_CV1se = fit_loto_glm$omega
)

for (method_name in names(omega_map)) {
  saveRDS(
    omega_map[[method_name]],
    file.path(omega_dir, sprintf("%s.rds", gsub("[^A-Za-z0-9]+", "_", method_name)))
  )
}

edge_dt <- rbindlist(lapply(names(omega_map), function(method_name) {
  omega_to_edge_dt(omega_map[[method_name]], taxa = colnames(X_raw), method = method_name)
}), fill = TRUE)

node_dt <- data.table(
  taxon = colnames(X_raw),
  total_abundance = colSums(X_raw),
  prevalence = colMeans(X_raw > 0),
  mean_log1p = colMeans(X_log),
  abundance_rank = seq_len(ncol(X_raw))
)

summary_dt <- data.table(
  method = names(omega_map),
  n_edges = vapply(omega_map, function(omega) sum(omega[upper.tri(omega)] != 0), numeric(1))
)

fwrite(edge_dt, file.path(out_dir, "network_full_fit_edges.csv"))
fwrite(node_dt, file.path(out_dir, "network_full_fit_nodes.csv"))
fwrite(summary_dt, file.path(out_dir, "network_full_fit_summary.csv"))

cat("Saved:\n")
cat(file.path(out_dir, "network_full_fit_edges.csv"), "\n")
cat(file.path(out_dir, "network_full_fit_nodes.csv"), "\n")
cat(file.path(out_dir, "network_full_fit_summary.csv"), "\n")
