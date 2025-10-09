library(data.table)
library(mlr3)
library(mlr3learners)
library(mlr3misc)
library(paradox)
library(checkmate)
library(glmnet)
library(ggplot2)
library(MASS)
library(broom)
library(statmod)
library(mlr3resampling)
library(torch)
library(batchtools)
library(igraph)
library(rlang)
library(PLNmodels)


dataname <- "glne007"
source("/projects/genomic-ml/da2343/PLN/pln_eval/load_source.R")
#task.dt <- data.table::fread("~/Projects/pln_eval/data/HMPv13_filtered.csv")
#task.dt <- data.table::fread(paste0("/projects/genomic-ml/da2343/PLN/pln_eval/data/", dataname, "_filtered.csv"))
task.dt <- data.table::fread(paste0("/projects/genomic-ml/da2343/PLN/pln_eval/data/", dataname, "_update.csv"))
#n_samples <- 30, 50, 100, 200, 500, 1000
n_samples <- 490
task.dt <- task.dt[1:n_samples,]
taxa_columns <- setdiff(names(task.dt), "Group_ID")
task.dt[, (taxa_columns) := lapply(.SD, function(x) log1p(x)), .SDcols = taxa_columns]

new_column_names <- paste0("Taxa", taxa_columns)
data.table::setnames(task.dt, old = taxa_columns, new = new_column_names)
task.list <- list()
for (col_name in new_column_names) {
  task_id <- col_name
  reg.task <- mlr3::TaskRegr$new(
    task_id, task.dt, target=col_name
  )
  reg.task$col_roles$feature <- setdiff(names(task.dt), c(task_id, "Group_ID"))
  task.list[[task_id]] <- reg.task
}


cv_glmnet_learner <- mlr3learners::LearnerRegrCVGlmnet$new()
cv_glmnet_learner$param_set$values$alpha <- 1
cv_glmnet_learner$param_set$values$type.measure <- "deviance"
cv_glmnet_learner$param_set$values$family <- "poisson"

plnpca_learner <- LearnerRegrPLNPCA$new()
plnpca_learner$param_set$values$rho <- paradox::to_tune(levels = c(0.02, 0.05, 0.1, 0.2, 0.5, 0.7))
plnpca_learner$encapsulate(method = "evaluate", fallback = mlr3::LearnerRegrFeatureless$new())
plnpca_learner_tuned <- mlr3tuning::auto_tuner(
  tuner = mlr3tuning::TunerBatchGridSearch$new(),
  learner = plnpca_learner,
  resampling = mlr3::ResamplingInsample$new(),
  measure = poisson_measure
)


pln_learner <- LearnerRegrPLN$new()
pln_learner$encapsulate(method = "evaluate", fallback = mlr3::LearnerRegrFeatureless$new())

pln_network_learner <- LearnerRegrPLNnetwork$new()
pln_network_learner$param_set$values$rho <- paradox::to_tune(levels = c(0.02, 0.05, 0.1, 0.2, 0.5, 0.7))
pln_network_learner$encapsulate(method = "evaluate", fallback = mlr3::LearnerRegrFeatureless$new())
pln_network_learner_tuned <- mlr3tuning::auto_tuner(
  tuner = mlr3tuning::TunerBatchGridSearch$new(),
  learner = pln_network_learner,
  resampling = mlr3::ResamplingInsample$new(),
  measure = poisson_measure
)



if(F){
#glmnet_learner_gaussian <- mlr3learners::LearnerRegrCVGlmnet$new()
#glmnet_learner_gaussian$param_set$values$alpha <- 1
#glmnet_learner_gaussian$param_set$values$type.measure <- "mse"
#glmnet_learner_gaussian$param_set$values$family <- "gaussian"
#glmnet_learner_gaussian$id <- "regr.glmnet_gaussian"
  
grid_tuner <- mlr3tuning::TunerBatchGridSearch$new()
subtrain_valid_cv <- mlr3::ResamplingCV$new()
subtrain_valid_cv$param_set$values$folds <- 2
plnpca_learner$param_set$values$rho <- paradox::to_tune(levels = c(0.01, 0.1, 0.5, 0.7, 0.9))
  
plnpca_learner_tuned <- mlr3tuning::auto_tuner(
    tuner = grid_tuner,
    learner = plnpca_learner,
    resampling = subtrain_valid_cv,
    measure = poisson_measure
  )
  
pln_network_learner$param_set$values$rho <- paradox::to_tune(levels = c(0.1, 0.5))
pln_network_learner_tuned <- mlr3tuning::auto_tuner(
  tuner = grid_tuner,
  learner = pln_network_learner,
  resampling = subtrain_valid_cv,
  measure = poisson_measure
)
}

reg.learner.list <- list(
  mlr3::LearnerRegrFeatureless$new(),
  cv_glmnet_learner,
  pln_learner,
  plnpca_learner_tuned,
  pln_network_learner_tuned
  #pln_network_learner,
  #plnpca_learner,
  #pln_learner,
  #pln_network_learner,
)


if(F){
### For debugging
debug_cv <- mlr3::ResamplingCV$new()
debug_cv$param_set$values$folds <- 3
debug.grid <- mlr3::benchmark_grid(
  task.list["Taxa799024"],
  reg.learner.list,
  debug_cv
)
debug.result <- mlr3::benchmark(debug.grid)
debug.score.dt <- debug.result$score(poisson_measure)
#debug.score.dt <- debug.result$score(mlr3::msr("regr.rmse"))

aggregate_results <- debug.score.dt[, .(
  mean_deviance =mean(regr.poisson_deviance),
  sd_deviance = sd(regr.poisson_deviance ),
  n_iterations = .N
), by = .(learner_id)]


print(aggregate_results)
  
# Set up cross-validation
mycv <- mlr3::ResamplingCV$new()
mycv$param_set$values$folds <- 3
reg.bench.grid <- mlr3::benchmark_grid(
  tasks = task.list,
  learners = reg.learner.list,
  resamplings = mycv
)
bmr <- mlr3::benchmark(reg.bench.grid)
results <- bmr$score(poisson_measure)
# Save results
save(bmr, file= paste0("~/Projects/pln_eval/data/", dataname, "_corr_benchmark.RData") )
# Extract only the simple columns for CSV export
results_clean <- results[, .(
  task_id, 
  learner_id, 
  iteration,
  regr.poisson_deviance
)]
write.csv(results_clean, paste0("~/Projects/pln_eval/data/", dataname, "_corr_benchmark.csv"), row.names = FALSE)
print(aggregate_results)
}

### For parallel real tests
future::plan("sequential")
mycv <- mlr3::ResamplingCV$new()
mycv$param_set$values$folds=3
(reg.bench.grid <- mlr3::benchmark_grid(
  task.list,
  reg.learner.list,
  mycv))

reg.dir <- paste0("/projects/genomic-ml/da2343/PLN/pln_eval/bmr/", dataname, "_", n_samples, "_samples")
#dir.create(reg.dir, recursive = TRUE)
reg <- batchtools::loadRegistry(reg.dir, writeable = TRUE)

#unlink(reg.dir, recursive=TRUE)

reg = batchtools::makeExperimentRegistry(
  file.dir = reg.dir,
  seed = 1,
  packages = c("mlr3verse", "batchtools", "data.table", "PLNmodels", "torch")
)
mlr3batchmark::batchmark(
  reg.bench.grid, store_models = FALSE, reg=reg)
job.table <- batchtools::getJobTable(reg=reg)
chunks <- data.frame(job.table, chunk=1)
batchtools::submitJobs(chunks, resources=list(
  walltime = 60*60*24,#seconds
  memory = 3072,#megabytes per cpu
  ncpus=1,  #>1 for multicore/parallel jobs.
  ntasks=1, #>1 for MPI jobs.
  nodelist = "cn67,cn65,cn63,cn62,cn61,cn58,cn56,cn55,cn54,cn52,cn51,cn50,cn49,cn48,cn47,cn46,cn44,cn41,cn36,cn23,cn22,cn21,cn19,cn17,cn16,cn15,cn14,cn13,cn12,cn4,cn101,cn99,cn98,cn95,cn93,cn92,cn90,cn89,cn88,cn87,cn83,cn82,cn81,cn80,cn79,cn78,cn77,cn71,cn106,cn30,cn29",
  chunks.as.arrayjobs=T), reg=reg)


batchtools::getStatus(reg=reg)
#batchtools::killJobs(reg=reg)
jobs.after <- batchtools::getJobTable(reg=reg)
table(jobs.after$error)
jobs.after[!is.na(error), .(error, task_id=sapply(prob.pars, "[[", "task_id"))][25:26]
#ids <- jobs.after[is.na(error), job.id]
ids <- jobs.after[!is.na(done) & is.na(error), job.id]
bmr = mlr3batchmark::reduceResultsBatchmark(ids, reg = reg)
score.dt <- bmr$score(poisson_measure)
#score.dt <- bmr$score(mlr3::msr("regr.rmse"))

aggregate_results <- score.dt[, .(
  mean_deviance = mean( regr.poisson_deviance  , na.rm = TRUE),
  sd_deviance = sd( regr.poisson_deviance , na.rm = TRUE),
  n_iterations = .N
), by = .(learner_id)]
aggregate_results[, n_samples := n_samples]

fwrite(aggregate_results, paste0("/projects/genomic-ml/da2343/PLN/pln_eval/bmr/", dataname, "_", n_samples, "_samples.csv"))
save(bmr, file=paste0("/projects/genomic-ml/da2343/PLN/pln_eval/bmr/", dataname, "_", n_samples, "_samples.RData"))


###############################
source("/projects/genomic-ml/da2343/PLN/pln_eval/load_source.R")
load("/projects/genomic-ml/da2343/PLN/pln_eval/bmr/amgut2_297_samples.RData")
score.dt <- as.data.table(bmr$score(poisson_measure))
# Map learner IDs to clean names
score.dt[, algorithm := fcase(
  learner_id == "regr.cv_glmnet", "GLMNet",
  learner_id == "regr.pln",       "PLN",
  default = NA_character_
)]
score.dt <- score.dt[!is.na(algorithm)]
score.dt <- score.dt[, .(
  algorithm,
  nr,
  iteration,
  deviance = regr.poisson_deviance
)]
score.dt[, n_samples := 297]
fwrite(score.dt, "/projects/genomic-ml/da2343/PLN/pln_eval/bmr/amgut2_297_09_11.csv")


