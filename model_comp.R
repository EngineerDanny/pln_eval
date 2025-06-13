library(data.table)
library(mlr3)
library(mlr3learners)
library(mlr3misc)
library(R6)
library(paradox)
library(checkmate)
library(glmnet)
library(ggplot2)

task.dt <- data.table::fread("/projects/genomic-ml/da2343/PLN/pln_eval/data/HMPv13_filtered.csv")
taxa_columns <- setdiff(names(task.dt), "Group_ID")
new_column_names <- paste0("Taxa", taxa_columns)
setnames(task.dt, old = taxa_columns, new = new_column_names)
task.list <- list()
for (col_name in new_column_names) {
  task_id <- col_name
  # task.list[[task_id]] <- mlr3::TaskRegr$new(
  # task_id, task.dt, target = col_name
  # )$set_col_roles("Group_ID", c("subset", "stratum"))
  
  reg.task <- mlr3::TaskRegr$new(
    task_id, task.dt, target=col_name
  )
  
  reg.task$col_roles$subset <- "Group_ID"
  # reg.task$col_roles$group <- "Group_ID"
  reg.task$col_roles$stratum <- "Group_ID"
  reg.task$col_roles$feature <- setdiff(names(task.dt), c(task_id, "Group_ID"))
  
  #task.list[[task_id]] <- mlr3::TaskRegr$new(
  #  task_id, task.dt, target=task_id
  #)$set_col_roles("Samples",c("subset", "stratum"))
  task.list[[task_id]] <- reg.task
}



set.seed(42)

# Then set up your learners with fallback
glmnet_learner <- mlr3learners::LearnerRegrCVGlmnet$new()
glmnet_learner$param_set$values$alpha <- 1
glmnet_learner$fallback <- mlr3::LearnerRegrFeatureless$new()
glmnet_learner$encapsulate <- c(train = "evaluate", predict = "evaluate")

reg.learner.list <- list(
  glmnet_learner,
  mlr3::LearnerRegrFeatureless$new()
)

## For debugging
debug_cv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
debug_cv$param_set$values$folds=5
(debug.grid <- mlr3::benchmark_grid(
  #task.list["TaxaTomentella"],
  task.list["Taxa4365684"],
  reg.learner.list,
  debug_cv))
debug.result <- mlr3::benchmark(debug.grid)
debug.score.dt <- mlr3resampling::score(debug.result)
aggregate_results <- debug.score.dt[, .(
  mean_mse = mean(regr.mse, na.rm = TRUE),
  sd_mse = sd(regr.mse, na.rm = TRUE),
  n_iterations = .N
), by = .(learner_id, train.subsets)]

print(aggregate_results)
