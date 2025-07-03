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
library(PLNmodels)
library(torch)


source("~/Projects/pln_eval/load_source.R")
task.dt <- data.table::fread("~/Projects/pln_eval/data/HMPv13_filtered.csv")
task.dt <- task.dt[1:50, 1:150]
taxa_columns <- setdiff(names(task.dt), "Group_ID")

# Apply log1p transformation to all abundance columns
task.dt[, (taxa_columns) := lapply(.SD, function(x) log1p(x)), .SDcols = taxa_columns]
new_column_names <- paste0("Taxa", taxa_columns)
setnames(task.dt, old = taxa_columns, new = new_column_names)

task.list <- list()
for (col_name in new_column_names) {
  task_id <- col_name
  reg.task <- mlr3::TaskRegr$new(
    task_id, task.dt, target=col_name
  )
  reg.task$col_roles$feature <- setdiff(names(task.dt), c(task_id, "Group_ID"))
  task.list[[task_id]] <- reg.task
}


# Create and use the measure
poisson_measure <- MeasurePoissonDeviance$new()
mlr3::mlr_measures$add("regr.poisson_deviance", MeasurePoissonDeviance)

cv_glmnet_learner <- mlr3learners::LearnerRegrCVGlmnet$new()
cv_glmnet_learner$param_set$values$alpha <- 1
cv_glmnet_learner$param_set$values$type.measure <- "deviance"
cv_glmnet_learner$param_set$values$family <- "poisson"

plnpca_learner <- LearnerRegrPLNPCA$new()
#plnpca_learner$param_set$rho <- 0.01
reg.learner.list <- list(
  mlr3::LearnerRegrFeatureless$new(),
  cv_glmnet_learner,
  LearnerRegrPLN$new(),
  #LearnerRegrCVPLN$new()
  plnpca_learner
)

## For debugging
debug_cv <- mlr3::ResamplingCV$new()
debug_cv$param_set$values$folds <- 3
debug.grid <- mlr3::benchmark_grid(
  task.list["Taxa799024"],
  reg.learner.list,
  debug_cv
)
debug.result <- mlr3::benchmark(debug.grid)

debug.score.dt <- debug.result$score(poisson_measure)
aggregate_results <- debug.score.dt[, .(
  mean_deviance =mean(regr.poisson_deviance)  ,
  sd_deviance = sd(regr.poisson_deviance ) ,
  n_iterations = .N
), by = .(learner_id)]

print(aggregate_results)

if(F){
learner_id mean_deviance sd_deviance n_iterations
<char>         <num>       <num>        <int>
  1: regr.featureless     2.3563887  0.52468928            3
2:   regr.cv_glmnet     0.9304769  0.25994747            3
3:         regr.pln     0.4355760  0.08876514            3
4:      regr.plnpca     0.4928400  0.15957130            3
}

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
