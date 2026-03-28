library(data.table)

args <- commandArgs(trailingOnly = TRUE)
base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
dataset_spec <- if (length(args) >= 1) args[1] else "amgut2_update.csv"

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

dataset_to_tag <- function(dataset_arg) {
  if (grepl("\\.tsv\\.gz$", dataset_arg)) {
    gsub("\\.tsv\\.gz$", "", basename(dataset_arg))
  } else if (grepl("\\.csv$", dataset_arg)) {
    gsub("\\.csv$", "", basename(dataset_arg))
  } else {
    dataset_arg
  }
}

expected_methods <- c(
  "Baseline_diagonal",
  "PLN_glasso",
  "LOTO_PLN_glasso",
  "LOTO_glmnet_CV1se"
)
expected_reps <- 1:3
expected_grid <- CJ(rep = expected_reps, method = expected_methods)

collect_one_dataset <- function(dataset_arg) {
  dataset_tag <- dataset_to_tag(dataset_arg)
  out_dir <- file.path(base_dir, "out", "network_graph", dataset_tag)
  task_dir <- file.path(out_dir, "cv3_tasks")
  if (!dir.exists(task_dir)) stop("Task directory not found: ", task_dir)

  result_files <- list.files(task_dir, pattern = "^result_fold.*\\.csv$", full.names = TRUE)
  edge_files <- list.files(task_dir, pattern = "^edges_fold.*\\.csv$", full.names = TRUE)
  if (length(result_files) == 0L) stop("No result files found in: ", task_dir)

  all_results <- rbindlist(lapply(result_files, fread), fill = TRUE)
  edge_dt <- if (length(edge_files) > 0L) {
    rbindlist(lapply(edge_files, fread), fill = TRUE)
  } else {
    data.table()
  }

  missing_dt <- expected_grid[!all_results, on = .(rep, method)]
  if (nrow(missing_dt) > 0L) {
    stop(
      "Missing CV task outputs for ", dataset_tag, ":\n",
      paste(sprintf("fold %d | %s", missing_dt$rep, missing_dt$method), collapse = "\n")
    )
  }

  summary_dt <- all_results[, .(
    mean_loglik = round(mean(loglik), 3),
    sd_loglik = round(sd(loglik), 3),
    mean_edges = round(mean(n_edges), 1)
  ), by = method][order(-mean_loglik)]

  fwrite(all_results, file.path(out_dir, "network_graph_pseudologlik_replicates.csv"))
  fwrite(summary_dt, file.path(out_dir, "network_graph_pseudologlik_summary.csv"))
  fwrite(edge_dt, file.path(out_dir, "network_graph_edges.csv"))

  cat("Saved for ", dataset_tag, ":\n", sep = "")
  cat(file.path(out_dir, "network_graph_pseudologlik_replicates.csv"), "\n")
  cat(file.path(out_dir, "network_graph_pseudologlik_summary.csv"), "\n")
  cat(file.path(out_dir, "network_graph_edges.csv"), "\n")
}

dataset_values <- read_dataset_spec(dataset_spec)
invisible(lapply(dataset_values, collect_one_dataset))
