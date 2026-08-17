#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(SpiecEasi)
})

args <- commandArgs(trailingOnly = TRUE)
datasets <- if (length(args)) args else c(
  "omm12",
  "omm12_keystone_2023",
  "pairinterax",
  "butyrate_assembly_2021",
  "host_fitness_2018"
)

base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
truth_dir <- file.path(base_dir, "data", "interaction_ground_truth")
out_path <- Sys.getenv(
  "SPIECEASI_OUTPUT",
  file.path(base_dir, "out", "spieceasi_bootstrap_f1_results.csv")
)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
B <- as.integer(Sys.getenv("SPIECEASI_BOOTSTRAPS", "20"))

read_tsv_gz <- function(path) {
  read.delim(
    gzfile(path), sep = "\t", header = TRUE, stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

sort_pair <- function(a, b) {
  ifelse(a <= b, paste(a, b, sep = "||"), paste(b, a, sep = "||"))
}

prepare_count_like <- function(X, target_depth = 10000L) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  X[!is.finite(X) | X < 0] <- 0
  if (all(abs(X - round(X)) < 1e-8)) return(round(X))

  totals <- rowSums(X)
  scaled <- matrix(0, nrow(X), ncol(X), dimnames = dimnames(X))
  keep <- totals > 0
  scaled[keep, ] <- round(X[keep, , drop = FALSE] / totals[keep] * target_depth)
  scaled
}

load_dataset <- function(dataset) {
  proc_dir <- file.path(truth_dir, dataset, "processed")
  abundance_file <- if (dataset %in% c("omm12", "omm12_keystone_2023")) {
    "community_absabundance_in_vitro.tsv.gz"
  } else {
    "abundance_matrix.tsv.gz"
  }
  abundance <- read_tsv_gz(file.path(proc_dir, abundance_file))
  X <- as.matrix(abundance[, -1, drop = FALSE])
  rownames(X) <- abundance$sample_id
  truth <- read_tsv_gz(file.path(proc_dir, "truth_undirected.tsv.gz"))
  list(X = X, truth = truth)
}

edge_f1 <- function(truth, adjacency, taxa) {
  truth_keys <- unique(sort_pair(truth$taxon_1, truth$taxon_2))
  edge_idx <- which(upper.tri(adjacency) & adjacency != 0, arr.ind = TRUE)
  pred_keys <- if (nrow(edge_idx)) {
    unique(sort_pair(taxa[edge_idx[, 1]], taxa[edge_idx[, 2]]))
  } else {
    character()
  }
  tp <- length(intersect(truth_keys, pred_keys))
  fp <- length(setdiff(pred_keys, truth_keys))
  fn <- length(setdiff(truth_keys, pred_keys))
  precision <- if (tp + fp > 0) tp / (tp + fp) else 0
  recall <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
  f1 <- if (tp == 0) {
    0
  } else if (is.finite(precision + recall) && precision + recall > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }
  list(f1 = f1, n_pred = length(pred_keys), tp = tp, fp = fp, fn = fn)
}

fit_spieceasi <- function(X) {
  fit <- spiec.easi(
    X,
    method = "glasso",
    nlambda = 50,
    lambda.min.ratio = 0.05,
    verbose = identical(Sys.getenv("SPIECEASI_VERBOSE"), "1")
  )
  opt_index <- getOptInd(fit)
  lambda_path <- fit$lambda
  list(
    fit = fit,
    adjacency = as.matrix(getRefit(fit)),
    selected_lambda = getOptLambda(fit),
    min_lambda = min(lambda_path),
    selected_min_lambda = opt_index == length(lambda_path)
  )
}

existing <- if (file.exists(out_path)) fread(out_path) else data.table()
rows <- list()

for (dataset_idx in seq_along(datasets)) {
  current_dataset <- datasets[dataset_idx]
  dat <- load_dataset(current_dataset)
  X <- prepare_count_like(dat$X)
  X <- X[rowSums(X) > 0, colSums(X) > 0, drop = FALSE]
  taxa <- colnames(X)

  for (b in seq_len(B)) {
    if (nrow(existing) && any(
      existing$dataset == current_dataset & existing$bootstrap_rep == b
    )) next
    seed <- 420000L + match(current_dataset, c(
      "omm12", "omm12_keystone_2023", "pairinterax",
      "butyrate_assembly_2021", "host_fitness_2018"
    )) * 1000L + b
    set.seed(seed)
    idx <- sample.int(nrow(X), nrow(X), replace = TRUE)
    X_boot <- X[idx, , drop = FALSE]
    rownames(X_boot) <- paste0("bs_", seq_len(nrow(X_boot)))

    started <- proc.time()[["elapsed"]]
    result <- tryCatch({
      fitted <- fit_spieceasi(X_boot)
      metrics <- edge_f1(dat$truth, fitted$adjacency, taxa)
      data.table(
        dataset = current_dataset, method = "SPIEC-EASI", bootstrap_rep = b,
        f1 = metrics$f1, n_pred = metrics$n_pred, tp = metrics$tp,
        fp = metrics$fp, fn = metrics$fn,
        selected_lambda = fitted$selected_lambda,
        min_lambda = fitted$min_lambda,
        selected_min_lambda = fitted$selected_min_lambda,
        elapsed_seconds = proc.time()[["elapsed"]] - started,
        seed = seed, error = NA_character_
      )
    }, error = function(e) {
      data.table(
        dataset = current_dataset, method = "SPIEC-EASI", bootstrap_rep = b,
        f1 = NA_real_, n_pred = NA_integer_, tp = NA_integer_,
        fp = NA_integer_, fn = NA_integer_,
        selected_lambda = NA_real_, min_lambda = NA_real_,
        selected_min_lambda = NA,
        elapsed_seconds = proc.time()[["elapsed"]] - started,
        seed = seed, error = conditionMessage(e)
      )
    })
    rows[[length(rows) + 1L]] <- result
    existing <- rbindlist(list(existing, result), fill = TRUE)
    setorder(existing, dataset, bootstrap_rep)
    fwrite(existing, out_path)
    message(
      current_dataset, " bootstrap ", b, "/", B, ": F1=",
      format(result$f1, digits = 3), ", edges=", result$n_pred,
      ", seconds=", round(result$elapsed_seconds, 1)
    )
  }
}

print(existing[, .(
  completed = sum(is.finite(f1)),
  mean_f1 = mean(f1, na.rm = TRUE),
  se_f1 = sd(f1, na.rm = TRUE) / sqrt(sum(is.finite(f1))),
  mean_seconds = mean(elapsed_seconds, na.rm = TRUE)
), by = dataset])
