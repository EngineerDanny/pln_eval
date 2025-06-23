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


task.dt <- data.table::fread("~/Projects/pln_eval/data/HMPv13_filtered.csv")
#task.dt <- task.dt[1:100]
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


LearnerRegrPLN <- R6::R6Class("LearnerRegrPLN",
                          inherit = LearnerRegr,
                          public = list(
                            initialize = function() {
                              ps = ps(
                                covariance = p_fct(c("full", "diagonal", "spherical"), default = "full", tags = "train"),
                                trace = p_int(0, 2, default = 0, tags = "train"),
                                backend = p_fct(c("nlopt", "torch"), default = "nlopt", tags = "train"),
                                offset_scheme = p_fct(c("TSS", "CSS", "RLE", "GMPR", "none"), default = "TSS", tags = "train")
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
                              pln_data <- prepare_data(
                                counts = abundance_matrix, 
                                covariates = covariates,
                                offset = pv$offset_scheme %||% "TSS"
                              )
                              
                              # Set up PLN control parameters
                              control_params <- PLN_param(
                                covariance = pv$covariance %||% "full",
                                trace = pv$trace %||% 0,
                                backend = "torch"
                              )
                              

                              
                              # Fit PLN model on prepared data
                              pln_model <- PLN(Abundance ~ 1, data = pln_data, control = control_params)
                              
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
                                test_pln_data <- prepare_data(
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
                                pln_cond_predictions <- predict_cond(
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
  mlr3::LearnerRegrFeatureless$new(),
  glmnet_learner,
  LearnerRegrPLN$new()
)

## For debugging
debug_cv <- mlr3::ResamplingCV$new()
debug_cv$param_set$values$folds <- 5
debug.grid <- mlr3::benchmark_grid(
  task.list["Taxa799024"],
  reg.learner.list,
  debug_cv
)
debug.result <- mlr3::benchmark(debug.grid)

debug.score.dt <- debug.result$score(poisson_measure)
aggregate_results <- debug.score.dt[, .(
  mean_deviance =mean(regr.poisson_deviance)  ,
  sd_deviance = sd(  regr.poisson_deviance ) ,
  n_iterations = .N
), by = .(learner_id)]

print(aggregate_results)



# Set up cross-validation
mycv <- mlr3::ResamplingCV$new()
mycv$param_set$values$folds <- 5
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
