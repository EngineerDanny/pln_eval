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
if (length(args) < 1) stop("Usage: Rscript model_comp.R <dataname>")
dataname <- args[1]

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
slurm_walltime <- 60L * 60L * 24L  # seconds
slurm_memory   <- 16384L            # MB per cpu

load_counts <- function(dataname, data_dir) {
  tsv_path <- list.files(data_dir, pattern = paste0(dataname, ".*\\.tsv\\.gz$"),
                         recursive = TRUE, full.names = TRUE)
  csv_path <- file.path(data_dir, paste0(dataname, "_update.csv"))

  if (length(tsv_path) > 0) {
    # taxa x samples → transpose to samples x taxa
    raw <- read.table(gzfile(tsv_path[1]), sep = "\t", header = TRUE,
                      row.names = 1, check.names = FALSE)
    counts <- t(as.matrix(raw))
    counts <- counts[rowSums(counts) > 0, colSums(counts) > 0]
    dt <- as.data.table(counts)
  } else if (file.exists(csv_path)) {
    dt <- data.table::fread(csv_path)
    if (!"Group_ID" %in% names(dt)) stop("Missing required column: Group_ID")
    dt[, Group_ID := NULL]
  } else {
    stop("No data file found for: ", dataname)
  }
  dt
}

task.dt <- load_counts(dataname, data_dir)

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
job.table <- batchtools::getJobTable(reg = reg)
chunks <- data.frame(job.table, chunk = 1)
batchtools::submitJobs(chunks, resources = list(
  walltime = slurm_walltime,
  memory   = slurm_memory,
  ncpus    = 1L,
  ntasks   = 1L,
  chunks.as.arrayjobs = TRUE
), reg = reg)

batchtools::waitForJobs(reg = reg)
batchtools::getStatus(reg = reg)
jobs.after <- batchtools::getJobTable(reg = reg)
table(jobs.after$error)

ids <- jobs.after[!is.na(done) & is.na(error), job.id]
bmr <- mlr3batchmark::reduceResultsBatchmark(ids, reg = reg)
score.dt <- bmr$score(poisson_measure)

aggregate_results <- score.dt[, .(
  mean_deviance = mean(regr.poisson_deviance, na.rm = TRUE),
  sd_deviance   = sd(regr.poisson_deviance, na.rm = TRUE),
  n_iterations  = .N
), by = .(learner_id)]
fwrite(aggregate_results, file.path(bmr_dir, paste0(dataname, ".csv")))
save(bmr, file = file.path(bmr_dir, paste0(dataname, ".RData")))
