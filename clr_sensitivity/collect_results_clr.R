library(data.table)
library(batchtools)
library(mlr3batchmark)
library(mlr3verse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript collect_results_clr.R <dataname>")

dataname <- args[1]
base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
registry_dir <- file.path(base_dir, "bmr_clr", dataname)
source(file.path(base_dir, "load_source.R"))
source(file.path(base_dir, "clr_sensitivity", "clr_learner.R"))

registry <- loadRegistry(registry_dir, writeable = FALSE)
waitForJobs(reg = registry)
jobs <- getJobTable(reg = registry)
if (any(!is.na(jobs$error))) print(jobs[!is.na(error), .(job.id, error)])

ids <- jobs[!is.na(done) & is.na(error), job.id]
if (!length(ids)) stop("No completed jobs found for ", dataname)

bmr <- reduceResultsBatchmark(ids, reg = registry)
scores <- as.data.table(bmr$score(poisson_measure))
summary <- scores[, .(
  mean_deviance = mean(regr.poisson_deviance, na.rm = TRUE),
  se_deviance = sd(regr.poisson_deviance, na.rm = TRUE) / sqrt(.N),
  n_iterations = .N
), by = learner_id]

fwrite(summary, file.path(base_dir, "bmr_clr", paste0(dataname, ".csv")))
save(bmr, file = file.path(base_dir, "bmr_clr", paste0(dataname, ".RData")))
print(summary)
