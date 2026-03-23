library(data.table)
library(mlr3)
library(mlr3misc)
library(paradox)
library(checkmate)
library(ggplot2)
library(mlr3resampling)
library(batchtools)
library(PLNmodels)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript model_comp_subsample.R <dataname_or_path> <n_samples> [memory_mb] [walltime_hours] [seed]")
dataname_or_path <- args[1]
n_sub            <- as.integer(args[2])
slurm_memory     <- if (length(args) >= 3) as.integer(args[3]) else 2048L
slurm_walltime   <- if (length(args) >= 4) as.integer(args[4]) * 3600L else 3600L
subsample_seed   <- if (length(args) >= 5) as.integer(args[5]) else 1L

# ----------------------------------------------------------------
# Paths
# ----------------------------------------------------------------
base_dir   <- "/projects/genomic-ml/da2343/PLN/pln_eval"
data_dir   <- file.path(base_dir, "data")
bmr_dir    <- file.path(base_dir, "bmr")

source(file.path(base_dir, "load_source.R"))

# ----------------------------------------------------------------
# Config
# ----------------------------------------------------------------
cv_folds <- 3L

load_counts <- function(dataname_or_path, data_dir) {
  if (grepl("\\.tsv\\.gz$", dataname_or_path)) {
    if (file.exists(dataname_or_path)) {
      tsv_path <- dataname_or_path
    } else {
      tsv_path <- file.path(data_dir, dataname_or_path)
    }
    if (!file.exists(tsv_path)) stop("File not found: ", dataname_or_path)
    raw <- read.table(gzfile(tsv_path), sep = "\t", header = TRUE,
                      row.names = 1, check.names = FALSE)
    counts <- t(as.matrix(raw))
    counts <- counts[rowSums(counts) > 0, colSums(counts) > 0]
    return(as.data.table(counts))
  }

  tsv_path <- list.files(data_dir, pattern = paste0(dataname_or_path, ".*\\.tsv\\.gz$"),
                         recursive = TRUE, full.names = TRUE)
  csv_path <- file.path(data_dir, paste0(dataname_or_path, "_update.csv"))

  if (length(tsv_path) > 0) {
    raw <- read.table(gzfile(tsv_path[1]), sep = "\t", header = TRUE,
                      row.names = 1, check.names = FALSE)
    counts <- t(as.matrix(raw))
    counts <- counts[rowSums(counts) > 0, colSums(counts) > 0]
    dt <- as.data.table(counts)
  } else if (file.exists(csv_path)) {
    dt <- data.table::fread(csv_path)
    if ("Group_ID" %in% names(dt)) dt[, Group_ID := NULL]
    dt
  } else {
    stop("No data file found for: ", dataname_or_path)
  }
}

# Derive clean name
dataname <- if (grepl("\\.tsv\\.gz$", dataname_or_path)) {
  gsub("\\.tsv\\.gz$", "", basename(dataname_or_path))
} else {
  dataname_or_path
}

task.dt <- load_counts(dataname_or_path, data_dir)
n_full <- nrow(task.dt)
p_full <- ncol(task.dt)

if (n_sub >= n_full) {
  stop(sprintf("Requested subsample %d >= full dataset size %d. Use model_comp.R instead.", n_sub, n_full))
}

# Subsample rows
set.seed(subsample_seed)
keep_rows <- sample.int(n_full, n_sub)
task.dt <- task.dt[keep_rows, ]

# Remove taxa with all zeros after subsampling
nonzero_cols <- colSums(task.dt) > 0
task.dt <- task.dt[, ..nonzero_cols]
p_after <- ncol(task.dt)

reg_name <- sprintf("%s_%d_samples", dataname, n_sub)
cat(sprintf("Loaded %s: %d/%d samples, %d/%d taxa (after filtering)\n",
            reg_name, n_sub, n_full, p_after, p_full))

# ----------------------------------------------------------------
# Build tasks
# ----------------------------------------------------------------
task.dt[, names(task.dt) := lapply(.SD, log1p)]

new_column_names <- paste0("Taxa", names(task.dt))
data.table::setnames(task.dt, names(task.dt), new_column_names)
task.list <- list()
for (col_name in new_column_names) {
  reg.task <- mlr3::TaskRegr$new(col_name, task.dt, target = col_name)
  task.list[[col_name]] <- reg.task
}

reg.learner.list <- list(
  mlr3::LearnerRegrFeatureless$new(),
  LearnerRegrLasso$new(),
  LearnerRegrPLN$new()
)

# ----------------------------------------------------------------
# Benchmark via batchtools
# ----------------------------------------------------------------
future::plan("sequential")
mycv <- mlr3::ResamplingCV$new()
mycv$param_set$values$folds <- cv_folds
reg.bench.grid <- mlr3::benchmark_grid(task.list, reg.learner.list, mycv)

reg.dir <- file.path(bmr_dir, reg_name)

if (dir.exists(reg.dir)) unlink(reg.dir, recursive = TRUE)
reg <- batchtools::makeExperimentRegistry(
  file.dir = reg.dir,
  seed = 1,
  packages = c("mlr3verse", "batchtools", "data.table", "PLNmodels", "torch"),
  source   = file.path(base_dir, "load_source.R")
)

mlr3batchmark::batchmark(reg.bench.grid, store_models = FALSE, reg = reg)

batchtools::submitJobs(
  resources = list(
    walltime = slurm_walltime, memory = slurm_memory,
    ncpus = 1L, ntasks = 1L, chunks.as.arrayjobs = TRUE
  ), reg = reg)

batchtools::getStatus(reg = reg)
cat("Jobs submitted. Run collect_results.R", reg_name, "to reduce results.\n")
