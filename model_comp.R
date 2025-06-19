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



Sys.setenv(LD_LIBRARY_PATH = "/packages/anaconda3/2024.02/lib")

# If TRUE, load them directly:
dyn.load("/packages/anaconda3/2024.02/lib/libicui18n.so.73")
dyn.load("/packages/anaconda3/2024.02/lib/libicuuc.so.73") 
dyn.load("/packages/anaconda3/2024.02/lib/libicudata.so.73")

library(PLNmodels)

task.dt <- data.table::fread("/projects/genomic-ml/da2343/PLN/pln_eval/data/HMPv13_filtered.csv")
task.dt <- task.dt[1:100]
taxa_columns <- setdiff(names(task.dt), "Group_ID")

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



LearnerRegrPLN <- R6::R6Class("LearnerRegrPLN",
                          inherit = LearnerRegr,
                          public = list(
                            initialize = function() {
                              ps = ps( 
                                covariance = paradox::p_fct(c("full", "diagonal", "spherical"), default = "full", tags = "train"),
                                trace = paradox::p_int(0, 2, default = 0, tags = "train"),
                                backend = paradox::p_fct(c("nlopt", "torch"), default = "nlopt", tags = "train"),
                                offset_scheme = paradox::p_fct(c("TSS", "CSS", "RLE", "GMPR", "none"), default = "TSS", tags = "train")
                              )
                              ps$values = list(covariance = "full", trace = 0, backend = "nlopt", offset_scheme = "TSS")
                              super$initialize(
                                id = "regr.pln",
                                param_set = ps,
                                feature_types = c("integer", "numeric"),
                                label = "Poisson Log-Normal Model",
                                packages = c("PLNmodels", "mlr3learners")
                              )
                            }
                          ),
                          private = list(
                            .train = function(task) {
                              pv <- self$param_set$get_values(tags = "train")
                              
                              # Get data from MLR3 task
                              X_train <- task$data(cols = task$feature_names)
                              y_train <- task$data(cols = task$target_names)[[1]]
                              
                              # Create abundance matrix: target + features
                              target_name <- task$target_names
                              feature_names <- task$feature_names
                              
                              # Combine target and features to create full abundance matrix
                              abundance_matrix <- cbind(y_train, X_train)
                              colnames(abundance_matrix) <- c(target_name, feature_names)
                              
                              # Set row names to match MLR3 task row IDs
                              abundance_matrix <- data.matrix(abundance_matrix, rownames.force = TRUE)
                              
                              # Create covariates (intercept only)
                              covariates <- data.frame(Intercept = rep(1, nrow(abundance_matrix)))
                              rownames(covariates) <- rownames(abundance_matrix)
                              
                              # Prepare data for PLN using prepare_data
                              # This may remove empty samples (rowSums == 0)
                              pln_data <- PLNmodels::prepare_data(
                                counts = abundance_matrix, 
                                covariates = covariates,
                                offset = pv$offset_scheme %||% "TSS"
                              )
                              
                              # Set up PLN control parameters
                              control_params <- PLNmodels::PLN_param(
                                covariance = pv$covariance %||% "full",
                                trace = pv$trace %||% 0,
                                backend = pv$backend %||% "torch"
                                #backend = pv$backend %||% "nlopt"
                              )
                              
                              # Fit PLN model on prepared data
                              pln_model <- PLNmodels::PLN(Abundance ~ 1, data = pln_data, control = control_params)
                              
                              # Store model and related info
                              self$model <- list(
                                pln_model = pln_model,
                                target_name = target_name,
                                feature_names = feature_names,
                                pv = pv
                              )
                              self$model
                            },
                            .predict = function(task) {
                              # Get test data from MLR3 task
                              X_test <- task$data(cols = task$feature_names)
                              y_test <- task$data(cols = task$target_names)[[1]]
                              
                              # Create test abundance matrix: target + features
                              target_name <- self$model$target_name
                              feature_names <- self$model$feature_names
                              
                              test_abundance_matrix <- cbind(y_test, X_test)
                              colnames(test_abundance_matrix) <- c(target_name, feature_names)
                              test_abundance_matrix <- data.matrix(test_abundance_matrix, rownames.force = TRUE)
                              
                              # Store original test info for MLR3 compatibility
                              n_test_original <- nrow(test_abundance_matrix)
                              original_row_names <- rownames(test_abundance_matrix)
                              
                              # Check for empty samples (rowSums == 0) before prepare_data
                              row_sums <- rowSums(test_abundance_matrix)
                              empty_samples <- which(row_sums == 0)
                              non_empty_samples <- which(row_sums > 0)
                              
                              # Initialize prediction vector
                              y_predict <- rep(NA_real_, n_test_original)
                              names(y_predict) <- original_row_names
                              
                              # Only process non-empty samples through PLN
                              if (length(non_empty_samples) > 0) {
                                # Create subset for non-empty samples
                                non_empty_abundance <- test_abundance_matrix[non_empty_samples, , drop = FALSE]
                                
                                # Create corresponding covariates
                                test_covariates <- data.frame(Intercept = rep(1, nrow(non_empty_abundance)))
                                rownames(test_covariates) <- rownames(non_empty_abundance)
                                
                                # Prepare data for PLN (should not remove any samples now)
                                test_pln_data <- PLNmodels::prepare_data(
                                  counts = non_empty_abundance,
                                  covariates = test_covariates,
                                  offset = self$model$pv$offset_scheme %||% "TSS"
                                )
                                
                                # For conditional prediction:
                                # newX = covariates (intercept)
                                # conditioning_data = feature species only (exclude target)
                                n_valid <- nrow(test_pln_data$Abundance)
                                newX <- data.frame("(Intercept)" = rep(1, n_valid), check.names = FALSE)
                                rownames(newX) <- rownames(test_pln_data$Abundance)
                                
                                # Conditioning data: feature species only (exclude target)
                                conditioning_data <- test_pln_data$Abundance[, feature_names, drop = FALSE]
                                
                                # Make conditional predictions using predict_cond
                                pln_cond_predictions <- PLNmodels::predict_cond(
                                  self$model$pln_model, 
                                  newX, 
                                  conditioning_data, 
                                  type = "response"
                                )
                                
                                # Extract prediction for target species (first column)
                                valid_predictions <- pln_cond_predictions[, 1]
                                
                                # Map predictions back to original indices
                                valid_row_names <- rownames(test_pln_data$Abundance)
                                y_predict[valid_row_names] <- valid_predictions
                              }
                              
                              # Handle empty samples with fallback
                              if (length(empty_samples) > 0) {
                                # For empty samples, predict zero counts
                                # (since they have no organisms, predicted count should be 0)
                                y_predict[empty_samples] <- 0
                              }
                              
                              # Ensure all predictions are filled (should not have NAs)
                              if (any(is.na(y_predict))) {
                                # Final fallback for any remaining NAs
                                remaining_nas <- is.na(y_predict)
                                fallback_value <- if (length(non_empty_samples) > 0) {
                                  mean(y_predict[!remaining_nas], na.rm = TRUE)
                                } else {
                                  0  # If all samples are empty, predict 0 counts
                                }
                                y_predict[remaining_nas] <- fallback_value
                              }
                              
                              # Return predictions (exactly matches original test sample count)
                              list(response = as.vector(y_predict))
                            }
                          )
)


MeasurePoissonDeviance <- R6::R6Class("MeasurePoissonDeviance",
                                      inherit = mlr3::MeasureRegr,
                                      public = list(
                                        initialize = function() {
                                          super$initialize(
                                            id = "regr.poisson_deviance",
                                            range = c(0, Inf),
                                            minimize = TRUE,
                                            predict_type = "response",
                                            packages = character(0),
                                            properties = character(0),
                                            label = "Poisson Deviance",
                                            man = "custom::poisson_deviance"
                                          )
                                        }
                                      ),
                                      private = list(
                                        .score = function(prediction, ...) {
                                          truth <- prediction$truth
                                          response <- prediction$response
                                          
                                          eps <- 1e-10
                                          response <- pmax(response, eps)
                                          log_term <- ifelse(truth == 0, 0, truth * log(truth / response))
                                          2 * mean(log_term - (truth - response))
                                        }
                                      )
)



# Create and use the measure
poisson_measure <- MeasurePoissonDeviance$new()
mlr3::mlr_measures$add("regr.poisson_deviance", MeasurePoissonDeviance)

glmnet_learner <- mlr3learners::LearnerRegrCVGlmnet$new()
glmnet_learner$param_set$values$alpha <- 1
glmnet_learner$param_set$values$type.measure <- "deviance"
glmnet_learner$param_set$values$family <- "poisson"
#glmnet_learner$fallback <- mlr3::LearnerRegrFeatureless$new()
#glmnet_learner$encapsulate <- c(train = "evaluate", predict = "evaluate")

reg.learner.list <- list(
  glmnet_learner,
  mlr3::LearnerRegrFeatureless$new(),
  LearnerRegrPLN$new()
)


### For debugging
debug_cv <- mlr3::ResamplingCV$new()
debug_cv$param_set$values$folds <- 5
debug.grid <- mlr3::benchmark_grid(
  task.list["Taxa4365684"],
  reg.learner.list,
  debug_cv
)
debug.result <- mlr3::benchmark(debug.grid)

debug.score.dt <- debug.result$score(poisson_measure)
#debug.score.dt <- debug.result$score(mlr3::msr("regr.rmse"))
aggregate_results <- debug.score.dt[, .(
  mean_deviance = mean( regr.poisson_deviance  , na.rm = TRUE),
  sd_deviance = sd( regr.poisson_deviance , na.rm = TRUE),
  n_iterations = .N
), by = .(learner_id, train.subsets)]
print(aggregate_results)



### For parallel real tests
future::plan("sequential")
mycv <- mlr3::ResamplingCV$new()
mycv$param_set$values$folds=5
(reg.bench.grid <- mlr3::benchmark_grid(
  task.list,
  reg.learner.list,
  mycv))

reg.dir <- "/projects/genomic-ml/da2343/PLN/pln_eval/out/hmpv13_06_19_reg"
reg <- batchtools::loadRegistry(reg.dir)
unlink(reg.dir, recursive=TRUE)
reg = batchtools::makeExperimentRegistry(
  file.dir = reg.dir,
  seed = 1,
  packages = "mlr3verse"
)
mlr3batchmark::batchmark(
  reg.bench.grid, store_models = TRUE, reg=reg)
job.table <- batchtools::getJobTable(reg=reg)
chunks <- data.frame(job.table, chunk=1)
batchtools::submitJobs(chunks, resources=list(
  walltime = 60*30,#seconds
  memory = 1024,#megabytes per cpu
  ncpus=1,  #>1 for multicore/parallel jobs.
  ntasks=1, #>1 for MPI jobs.
  chunks.as.arrayjobs=T), reg=reg)


batchtools::getStatus(reg=reg)
#batchtools::killJobs(reg=reg)
jobs.after <- batchtools::getJobTable(reg=reg)
table(jobs.after$error)
jobs.after[!is.na(error), .(error, task_id=sapply(prob.pars, "[[", "task_id"))][25:26]



ids <- jobs.after[is.na(error), job.id]
bmr = mlr3batchmark::reduceResultsBatchmark(ids, reg = reg)
score.dt <- mlr3resampling::score(bmr)
save(bmr, file="/projects/genomic-ml/da2343/PLN/pln_eval/out/hmpv13.RData")

jobs.final <- batchtools::getJobTable(reg=reg)
ids <- jobs.final[!is.na(done), job.id]
bmr = mlr3batchmark::reduceResultsBatchmark(ids, reg = reg)
score.dt <- mlr3resampling::score(bmr, poisson_measure)
aggregate_results <- score.dt[, .(
  mean_deviance = mean( regr.poisson_deviance  , na.rm = TRUE),
  sd_deviance = sd( regr.poisson_deviance , na.rm = TRUE),
  n_iterations = .N
), by = .(learner_id, train.subsets)]
print(aggregate_results)
