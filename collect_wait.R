library(data.table)
library(batchtools)
library(mlr3batchmark)
library(mlr3verse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript collect_wait.R <dataname>")
dataname <- args[1]

base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
bmr_dir  <- file.path(base_dir, "bmr")
reg.dir  <- file.path(bmr_dir, dataname)

source(file.path(base_dir, "load_source.R"))

reg <- batchtools::loadRegistry(reg.dir, writeable = FALSE)

cat("Waiting for jobs to finish for", dataname, "...\n")
batchtools::waitForJobs(reg = reg, sleep = 60)
cat("All jobs done.\n")

batchtools::getStatus(reg = reg)

jobs <- batchtools::getJobTable(reg = reg)
cat("Error table:\n")
print(table(jobs$error))

ids <- jobs[!is.na(done) & is.na(error), job.id]
cat("Completed jobs:", length(ids), "\n")

bmr      <- mlr3batchmark::reduceResultsBatchmark(ids, reg = reg)
score.dt <- bmr$score(poisson_measure)

aggregate_results <- score.dt[, .(
  mean_deviance = mean(regr.poisson_deviance, na.rm = TRUE),
  sd_deviance   = sd(regr.poisson_deviance, na.rm = TRUE),
  n_iterations  = .N
), by = .(learner_id)]

fwrite(aggregate_results, file.path(bmr_dir, paste0(dataname, ".csv")))
save(bmr, file = file.path(bmr_dir, paste0(dataname, ".RData")))

cat("\nResults for", dataname, ":\n")
print(aggregate_results)
