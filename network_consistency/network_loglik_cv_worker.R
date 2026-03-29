library(data.table)
library(PLNmodels)
library(glassoFast)
library(glmnet)

args <- commandArgs(trailingOnly = TRUE)
base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
data_dir <- file.path(base_dir, "data")

read_dataset_spec <- function(dataset_spec) {
  if (length(dataset_spec) != 1L || !nzchar(dataset_spec)) {
    stop("dataset_spec must be a non-empty string")
  }

  if (file.exists(dataset_spec) && grepl("\\.(txt|list|lst)$", dataset_spec)) {
    lines <- trimws(readLines(dataset_spec, warn = FALSE))
    lines <- lines[nzchar(lines)]
    lines <- lines[!grepl("^#", lines)]
    if (length(lines) == 0L) stop("No datasets found in: ", dataset_spec)
    return(lines)
  }

  dataset_spec
}

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

dataset_spec <- if (length(args) >= 1) args[1] else "amgut2_update.csv"
task_id <- if (length(args) >= 2) as.integer(args[2]) else {
  as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
}
if (!is.finite(task_id) || task_id < 1L) stop("Invalid task_id")

dataset_values <- read_dataset_spec(dataset_spec)
n_datasets <- length(dataset_values)
n_folds <- 3L
methods <- c(
  "Baseline_diagonal",
  "PLN_glasso",
  "LOTO_PLN_glasso",
  "LOTO_glmnet_CV1se"
)
n_methods <- length(methods)
n_tasks <- n_datasets * n_folds * n_methods
if (task_id > n_tasks) stop("task_id exceeds number of tasks: ", n_tasks)

dataset_id <- ((task_id - 1L) %/% (n_folds * n_methods)) + 1L
within_dataset_id <- ((task_id - 1L) %% (n_folds * n_methods)) + 1L
fold_id <- ((within_dataset_id - 1L) %/% n_methods) + 1L
method_id <- ((within_dataset_id - 1L) %% n_methods) + 1L

dataset_arg <- dataset_values[dataset_id]
method_name <- methods[method_id]

dataset_tag <- if (grepl("\\.tsv\\.gz$", dataset_arg)) {
  gsub("\\.tsv\\.gz$", "", basename(dataset_arg))
} else if (grepl("\\.csv$", dataset_arg)) {
  gsub("\\.csv$", "", basename(dataset_arg))
} else {
  dataset_arg
}

out_dir <- file.path(base_dir, "out", "network_graph", dataset_tag)
task_dir <- file.path(out_dir, "cv3_tasks")
dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)

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

omega_to_edge_dt <- function(Omega, taxa, fold_id, method) {
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
    rep = fold_id,
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

support_to_edge_dt <- function(Omega, support, taxa, fold_id, method) {
  idx <- which(upper.tri(support), arr.ind = TRUE)
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
    rep = fold_id,
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

empty_support <- function(p) {
  matrix(FALSE, p, p)
}

support_from_precision <- function(Omega, tol = 1e-8) {
  support <- abs(Omega) > tol
  diag(support) <- FALSE
  support | t(support)
}

support_from_coef_mat <- function(coef_mat, require_mutual = TRUE, tol = 1e-8) {
  nz <- abs(coef_mat) > tol
  diag(nz) <- FALSE
  support <- if (require_mutual) nz & t(nz) else nz | t(nz)
  diag(support) <- FALSE
  support
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
  list(omega = gl_best$wi, best_rho = best_rho)
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

evaluate_support_nodewise <- function(X_train, X_test, support, ridge_lambda = 1e-3) {
  p <- ncol(X_train)
  coef_mat <- matrix(0, p, p)
  improvement <- numeric(p)

  for (j in seq_len(p)) {
    nbrs <- which(support[j, ])
    y_train <- X_train[, j]
    y_test <- X_test[, j]
    baseline_pred <- mean(y_train)
    mse_base <- mean((y_test - baseline_pred)^2)

    if (length(nbrs) == 0L) {
      pred_test <- rep(baseline_pred, length(y_test))
    } else {
      Xj_train <- as.matrix(X_train[, nbrs, drop = FALSE])
      Xj_test <- as.matrix(X_test[, nbrs, drop = FALSE])
      x_mu <- colMeans(Xj_train)
      y_mu <- mean(y_train)
      Xj_train_centered <- sweep(Xj_train, 2, x_mu, "-")
      Xj_test_centered <- sweep(Xj_test, 2, x_mu, "-")
      y_train_centered <- y_train - y_mu

      XtX <- crossprod(Xj_train_centered)
      beta <- tryCatch(
        solve(
          XtX + diag(ridge_lambda, ncol(Xj_train)),
          crossprod(Xj_train_centered, y_train_centered)
        ),
        error = function(e) rep(0, ncol(Xj_train))
      )
      coef_mat[j, nbrs] <- as.numeric(beta)
      pred_test <- y_mu + as.vector(Xj_test_centered %*% beta)
    }

    mse_model <- mean((y_test - pred_test)^2)
    improvement[j] <- if (is.finite(mse_base) && mse_base > 1e-8) {
      1 - mse_model / mse_base
    } else {
      0
    }
  }

  list(
    score = mean(improvement),
    coef_mat = coef_mat
  )
}

counts_full <- load_counts(dataset_arg, data_dir)
col_sums <- colSums(counts_full)
top_taxa <- order(col_sums, decreasing = TRUE)
X_raw <- counts_full[, top_taxa, drop = FALSE]
X_log <- log1p(X_raw)
rownames(X_raw) <- paste0("S", seq_len(nrow(X_raw)))
rownames(X_log) <- rownames(X_raw)

n_total <- nrow(X_log)
p <- ncol(X_log)
cat(sprintf("Data: %s | %d samples, %d taxa (log1p)\n", dataset_tag, n_total, p))

set.seed(42)
fold_assign <- sample(rep(seq_len(n_folds), length.out = n_total))
train_idx <- which(fold_assign != fold_id)
test_idx <- which(fold_assign == fold_id)

X_train_raw <- X_raw[train_idx, , drop = FALSE]
X_train <- X_log[train_idx, , drop = FALSE]
X_test <- X_log[test_idx, , drop = FALSE]

cat(sprintf(
  "Task %d/%d | dataset=%d/%d | fold=%d | method=%s | train=%d | test=%d\n",
  task_id, n_tasks, dataset_id, n_datasets, fold_id, method_name, length(train_idx), length(test_idx)
))

t0 <- proc.time()
if (method_name == "Baseline_diagonal") {
  support <- empty_support(p)
} else if (method_name == "PLN_glasso") {
  fit <- fit_pln_glasso(X_train_raw)
  support <- support_from_precision(fit$omega)
} else if (method_name == "LOTO_PLN_glasso") {
  fit <- fit_loto_pln_glasso(X_train_raw, X_train)
  support <- support_from_coef_mat(fit$coef_mat, require_mutual = TRUE)
} else if (method_name == "LOTO_glmnet_CV1se") {
  coef_mat <- matrix(0, p, p)
  for (j in seq_len(p)) {
    y_train <- X_train_raw[, j]
    x_train <- X_train[, -j, drop = FALSE]
    coef_mat[j, -j] <- safe_glmnet_coefs(x_train, y_train)
  }
  support <- support_from_coef_mat(coef_mat, require_mutual = TRUE)
} else {
  stop("Unknown method: ", method_name)
}

eval_fit <- evaluate_support_nodewise(X_train, X_test, support)
n_edges <- sum(support[upper.tri(support)])
edge_weights <- (eval_fit$coef_mat + t(eval_fit$coef_mat)) / 2

score <- eval_fit$score
elapsed <- (proc.time() - t0)[["elapsed"]]

result_dt <- data.table(
  rep = fold_id,
  method = method_name,
  score = round(score, 6),
  n_edges = n_edges,
  elapsed_sec = elapsed
)
edge_dt <- support_to_edge_dt(edge_weights, support, colnames(X_train), fold_id, method_name)

result_path <- file.path(task_dir, sprintf("result_fold%02d_%s.csv", fold_id, method_name))
edge_path <- file.path(task_dir, sprintf("edges_fold%02d_%s.csv", fold_id, method_name))
fwrite(result_dt, result_path)
fwrite(edge_dt, edge_path)

cat("Saved:\n")
cat(result_path, "\n")
cat(edge_path, "\n")
