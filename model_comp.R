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
if (length(args) < 1) {
  stop("Usage: Rscript model_comp.R <dataname_or_path> [memory_mb] [walltime_hours]")
}
dataname_or_path <- args[1]
slurm_memory     <- if (length(args) >= 2) as.integer(args[2]) else 2048L    # MB per cpu
slurm_walltime   <- if (length(args) >= 3) as.integer(args[3]) * 3600L else 3600L  # default 1h
array_max_size   <- 25000L

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
rho_grid       <- c(0.02, 0.05, 0.1, 0.2, 0.5, 0.7)
cv_folds       <- 3L

load_counts <- function(dataname_or_path, data_dir) {
  # Direct path to tsv.gz (absolute or relative to data_dir)
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

  # Pattern match on dataname
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

# Derive a clean registry name from the input
dataname <- if (grepl("\\.tsv\\.gz$", dataname_or_path)) {
  gsub("\\.tsv\\.gz$", "", basename(dataname_or_path))
} else {
  dataname_or_path
}

task.dt <- load_counts(dataname_or_path, data_dir)
cat(sprintf("Loaded %s: %d samples x %d taxa\n", dataname, nrow(task.dt), ncol(task.dt)))

# ----------------------------------------------------------------
# Validation
# ----------------------------------------------------------------
if (nrow(task.dt) == 0) stop("Data is empty: no rows found in ", dataname)
if (ncol(task.dt) == 0) stop("Data is empty: no columns found in ", dataname)

task.dt[, names(task.dt) := lapply(.SD, log1p)]

new_column_names <- paste0("Taxa", names(task.dt))
data.table::setnames(task.dt, names(task.dt), new_column_names)
task.list <- list()
for (col_name in new_column_names) {
  reg.task <- mlr3::TaskRegr$new(col_name, task.dt, target = col_name)
  task.list[[col_name]] <- reg.task
}

lasso_learner <- LearnerRegrLasso$new()

# plnpca_learner <- LearnerRegrPLNPCA$new()
# plnpca_learner$param_set$values$rho <- paradox::to_tune(levels = rho_grid)
# plnpca_learner$encapsulate(method = "evaluate", fallback = mlr3::LearnerRegrFeatureless$new())
# plnpca_learner_tuned <- mlr3tuning::auto_tuner(
#   tuner = mlr3tuning::TunerBatchGridSearch$new(),
#   learner = plnpca_learner,
#   resampling = mlr3::ResamplingInsample$new(),
#   measure = poisson_measure
# )

pln_learner <- LearnerRegrPLN$new()

reg.learner.list <- list(
  mlr3::LearnerRegrFeatureless$new(),
  lasso_learner,
  pln_learner
)

### Parallel benchmark via batchtools
future::plan("sequential")
mycv <- mlr3::ResamplingCV$new()
mycv$param_set$values$folds <- cv_folds
reg.bench.grid <- mlr3::benchmark_grid(task.list, reg.learner.list, mycv)

reg.dir <- file.path(bmr_dir, dataname)

if (dir.exists(reg.dir)) unlink(reg.dir, recursive = TRUE)
reg <- batchtools::makeExperimentRegistry(
  file.dir = reg.dir,
  seed = 1,
  packages = c("mlr3verse", "batchtools", "data.table", "PLNmodels", "torch"),
  source   = file.path(base_dir, "load_source.R")
)

mlr3batchmark::batchmark(reg.bench.grid, store_models = FALSE, reg = reg)
job.table <- batchtools::findNotSubmitted(reg = reg)
job.table[, chunk := batchtools::chunk(job.id, chunk.size = array_max_size, shuffle = FALSE)]

cat(sprintf(
  "Preparing %d jobs in %d array submissions (max array size %d)\n",
  nrow(job.table), data.table::uniqueN(job.table$chunk), array_max_size
))

batchtools::submitJobs(
  ids = job.table,
  resources = list(
    walltime = slurm_walltime, memory = slurm_memory,
    ncpus = 1L, ntasks = 1L, chunks.as.arrayjobs = TRUE
  ), reg = reg)

batchtools::getStatus(reg = reg)
cat("Jobs submitted. Run collect_results.R", dataname, "to reduce results.\n")
