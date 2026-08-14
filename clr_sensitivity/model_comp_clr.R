library(data.table)
library(mlr3)
library(mlr3misc)
library(batchtools)
library(glmnet)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript model_comp_clr.R <dataname_or_path> [memory_mb] [walltime_hours]")
}

dataname_or_path <- args[1]
slurm_memory <- if (length(args) >= 2) as.integer(args[2]) else 2048L
slurm_walltime <- if (length(args) >= 3) as.integer(args[3]) * 3600L else 3600L
array_max_size <- 25000L

base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
data_dir <- file.path(base_dir, "data")
registry_root <- file.path(base_dir, "bmr_clr")

source(file.path(base_dir, "load_source.R"))
source(file.path(base_dir, "clr_sensitivity", "clr_learner.R"))

load_counts <- function(dataname_or_path, data_dir) {
  known_paths <- c(
    crc_zeller = "family/microbiomehd/crc_zeller_family.tsv.gz"
  )
  if (dataname_or_path %in% names(known_paths)) {
    dataname_or_path <- unname(known_paths[dataname_or_path])
  }

  if (grepl("\\.tsv\\.gz$", dataname_or_path)) {
    path <- if (file.exists(dataname_or_path)) {
      dataname_or_path
    } else {
      file.path(data_dir, dataname_or_path)
    }
    if (!file.exists(path)) stop("File not found: ", dataname_or_path)
    raw <- read.table(
      gzfile(path), sep = "\t", header = TRUE, row.names = 1,
      check.names = FALSE
    )
    counts <- t(as.matrix(raw))
    return(counts[rowSums(counts) > 0, colSums(counts) > 0, drop = FALSE])
  }

  matches <- list.files(
    data_dir,
    pattern = paste0("^", dataname_or_path, ".*\\.tsv\\.gz$"),
    recursive = TRUE,
    full.names = TRUE
  )
  csv_path <- file.path(data_dir, paste0(dataname_or_path, "_update.csv"))

  if (length(matches) == 1L) {
    raw <- read.table(
      gzfile(matches), sep = "\t", header = TRUE, row.names = 1,
      check.names = FALSE
    )
    counts <- t(as.matrix(raw))
    return(counts[rowSums(counts) > 0, colSums(counts) > 0, drop = FALSE])
  }
  if (length(matches) > 1L) {
    stop("Multiple data files matched ", dataname_or_path, "; pass a relative path instead")
  }
  if (file.exists(csv_path)) {
    counts <- fread(csv_path)
    if ("Group_ID" %in% names(counts)) counts[, Group_ID := NULL]
    return(as.matrix(counts))
  }
  stop("No data file found for: ", dataname_or_path)
}

dataname <- if (grepl("\\.tsv\\.gz$", dataname_or_path)) {
  sub("\\.tsv\\.gz$", "", basename(dataname_or_path))
} else {
  dataname_or_path
}

counts <- load_counts(dataname_or_path, data_dir)
storage.mode(counts) <- "double"
if (!nrow(counts) || !ncol(counts)) stop("Dataset is empty: ", dataname)

# Match the primary benchmark: log1p response and predictor abundances.
task_dt <- as.data.table(log1p(counts))
setnames(task_dt, names(task_dt), paste0("Taxa", names(task_dt)))
tasks <- lapply(names(task_dt), function(target) {
  TaskRegr$new(target, task_dt, target = target)
})
names(tasks) <- names(task_dt)

learners <- list(
  LearnerRegrFeatureless$new(),
  LearnerRegrLasso$new(),
  LearnerRegrCLRLasso$new()
)
resampling <- ResamplingCV$new()
resampling$param_set$values$folds <- 3L
grid <- benchmark_grid(tasks, learners, resampling)

dir.create(registry_root, recursive = TRUE, showWarnings = FALSE)
registry_dir <- file.path(registry_root, dataname)
if (dir.exists(registry_dir)) {
  stop("Registry already exists: ", registry_dir, ". Move or remove it before rerunning.")
}

future::plan("sequential")
registry <- makeExperimentRegistry(
  file.dir = registry_dir,
  seed = 1,
  packages = c("mlr3verse", "batchtools", "data.table", "glmnet"),
  source = c(
    file.path(base_dir, "load_source.R"),
    file.path(base_dir, "clr_sensitivity", "clr_learner.R")
  )
)
mlr3batchmark::batchmark(grid, store_models = FALSE, reg = registry)

jobs <- findNotSubmitted(reg = registry)
jobs[, chunk := chunk(job.id, chunk.size = array_max_size, shuffle = FALSE)]
submitJobs(
  ids = jobs,
  resources = list(
    walltime = slurm_walltime,
    memory = slurm_memory,
    ncpus = 1L,
    ntasks = 1L,
    chunks.as.arrayjobs = TRUE
  ),
  reg = registry
)

cat(sprintf(
  "Submitted %s: %d samples, %d taxa, %d jobs\n",
  dataname, nrow(task_dt), ncol(task_dt), nrow(jobs)
))
