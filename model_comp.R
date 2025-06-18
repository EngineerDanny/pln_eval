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



LearnerRegrPLNGlmnet <- R6::R6Class("LearnerRegrPLNGlmnet",
                                    inherit = LearnerRegr,
                                    
                                    public = list(
                                      
                                      #' @description Create a new combined PLN+glmnet learner
                                      initialize = function() {
                                        ps = ps(
                                          # cv.glmnet parameters for feature selection
                                          alpha = p_dbl(0, 1, default = 1, tags = "train"),
                                          lambda_choice = p_fct(c("lambda.min", "lambda.1se"), default = "lambda.min", tags = "train"),
                                          nfolds_glmnet = p_int(3, 10, default = 5, tags = "train"),
                                          family_glmnet = p_fct(c("poisson", "gaussian"), default = "poisson", tags = "train"),
                                          
                                          # PLN parameters for joint modeling
                                          covariance = p_fct(c("full", "diagonal", "spherical"), default = "full", tags = "train"),
                                          backend = p_fct(c("nlopt", "torch"), default = "torch", tags = "train"),
                                          trace = p_int(0, 2, default = 0, tags = "train"),
                                          offset_scheme = p_fct(c("TSS", "CSS", "RLE", "GMPR", "none"), default = "TSS", tags = "train"),
                                          
                                          # Combination strategy
                                          weight_glmnet = p_dbl(0, 1, default = 0.3, tags = "train"),  # Weight for glmnet predictions
                                          use_selected_features = p_lgl(default = TRUE, tags = "train")  # Whether to use feature selection info
                                        )
                                        
                                        ps$values = list(
                                          alpha = 1,
                                          lambda_choice = "lambda.min", 
                                          nfolds_glmnet = 5,
                                          family_glmnet = "poisson",
                                          covariance = "full",
                                          backend = "torch",
                                          trace = 0,
                                          offset_scheme = "TSS",
                                          weight_glmnet = 0.3,
                                          use_selected_features = TRUE
                                        )
                                        
                                        super$initialize(
                                          id = "regr.pln_glmnet_combined",
                                          param_set = ps,
                                          feature_types = c("integer", "numeric"),
                                          label = "Combined PLN + cv.glmnet Learner",
                                          packages = c("PLNmodels", "glmnet", "mlr3learners")
                                        )
                                      }
                                    ),
                                    
                                    private = list(
                                      
                                      #' @description Train the combined learner
                                      .train = function(task) {
                                        pv <- self$param_set$get_values(tags = "train")
                                        
                                        # Extract data from task
                                        X_features <- task$data(cols = task$feature_names)
                                        y_target <- task$data(cols = task$target_names)[[1]]
                                        target_name <- task$target_names
                                        feature_names <- task$feature_names
                                        
                                        # Create abundance matrix: target + features (as species)
                                        abundance_matrix <- cbind(y_target, X_features)
                                        colnames(abundance_matrix) <- c(target_name, feature_names)
                                        abundance_matrix <- as.matrix(abundance_matrix)
                                        rownames(abundance_matrix) <- paste0("sample_", seq_len(nrow(abundance_matrix)))
                                        
                                        # Create minimal covariates (intercept only for PLN)
                                        covariates_df <- data.frame(
                                          intercept = rep(1, nrow(abundance_matrix)),
                                          row.names = rownames(abundance_matrix)
                                        )
                                        
                                        # Prepare data using real prepare_data function
                                        pln_data <- prepare_data(
                                          counts = abundance_matrix,
                                          covariates = covariates_df,
                                          offset = pv$offset_scheme
                                        )
                                        
                                        ## =======================================================
                                        ## STEP 1: cv.glmnet for feature selection and baseline predictions
                                        ## =======================================================
                                        
                                        if (pv$trace > 0) message("Step 1: cv.glmnet feature selection...")
                                        
                                        glmnet_model <- NULL
                                        selected_features <- character(0)
                                        feature_importance <- NULL
                                        
                                        tryCatch({
                                          # Prepare data for glmnet: target ~ features  
                                          y_glmnet <- pln_data$Abundance[, target_name]
                                          X_glmnet <- pln_data$Abundance[, feature_names, drop = FALSE]
                                          
                                          # Fit cv.glmnet
                                          glmnet_model <- cv.glmnet(
                                            x = X_glmnet,
                                            y = y_glmnet, 
                                            family = pv$family_glmnet,
                                            alpha = pv$alpha,
                                            nfolds = pv$nfolds_glmnet
                                          )
                                          
                                          # Extract selected features and their importance
                                          coef_glmnet <- coef(glmnet_model, s = glmnet_model[[pv$lambda_choice]])
                                          nonzero_idx <- which(coef_glmnet[-1, 1] != 0)  # Exclude intercept
                                          
                                          if (length(nonzero_idx) > 0) {
                                            selected_features <- feature_names[nonzero_idx]
                                            feature_importance <- abs(coef_glmnet[-1, 1][nonzero_idx])
                                            names(feature_importance) <- selected_features
                                          }
                                          
                                          if (pv$trace > 0) {
                                            message(sprintf("cv.glmnet selected %d/%d features (%.1f%% sparsity)", 
                                                            length(selected_features), length(feature_names),
                                                            (1 - length(selected_features)/length(feature_names)) * 100))
                                          }
                                          
                                        }, error = function(e) {
                                          if (pv$trace > 0) message("cv.glmnet failed: ", e$message)
                                        })
                                        
                                        ## =======================================================
                                        ## STEP 2: PLN for joint modeling
                                        ## =======================================================
                                        
                                        if (pv$trace > 0) message("Step 2: PLN joint modeling...")
                                        
                                        pln_model <- NULL
                                        
                                        tryCatch({
                                          # Create PLN control parameters - use only real parameters
                                          control_pln <- list(
                                            covariance = pv$covariance,
                                            trace = pv$trace,
                                            backend = pv$backend
                                          )
                                          
                                          # Fit PLN model using real PLN function
                                          pln_model <- PLN(Abundance ~ 1, data = pln_data, control = control_pln)
                                          
                                          if (pv$trace > 0) {
                                            message(sprintf("PLN model fitted. Loglik: %.2f, BIC: %.2f", 
                                                            pln_model$loglik, pln_model$BIC))
                                          }
                                          
                                        }, error = function(e) {
                                          if (pv$trace > 0) message("PLN failed: ", e$message)
                                        })
                                        
                                        ## =======================================================
                                        ## STEP 3: Store combined model
                                        ## =======================================================
                                        
                                        # Calculate sparsity metrics
                                        sparsity_level <- ifelse(length(feature_names) > 0, 
                                                                 1 - length(selected_features) / length(feature_names), 
                                                                 0)
                                        
                                        self$model <- list(
                                          # Models
                                          glmnet_model = glmnet_model,
                                          pln_model = pln_model,
                                          
                                          # Data structure
                                          pln_data = pln_data,
                                          target_name = target_name,
                                          feature_names = feature_names,
                                          abundance_colnames = colnames(abundance_matrix),
                                          
                                          # Feature selection results
                                          selected_features = selected_features,
                                          feature_importance = feature_importance,
                                          sparsity_level = sparsity_level,
                                          
                                          # Parameters
                                          params = pv,
                                          
                                          # Model validation
                                          glmnet_success = !is.null(glmnet_model),
                                          pln_success = !is.null(pln_model)
                                        )
                                        
                                        if (pv$trace > 0) {
                                          message(sprintf("Combined model ready. GLMnet: %s, PLN: %s", 
                                                          ifelse(self$model$glmnet_success, "✓", "✗"),
                                                          ifelse(self$model$pln_success, "✓", "✗")))
                                        }
                                        
                                        return(self$model)
                                      },
                                      
                                      #' @description Make combined predictions
                                      .predict = function(task) {
                                        
                                        # Extract test data
                                        X_test <- task$data(cols = task$feature_names)
                                        y_test <- task$data(cols = task$target_names)[[1]]
                                        
                                        # Recreate abundance matrix structure
                                        test_abundance <- cbind(y_test, X_test)
                                        colnames(test_abundance) <- self$model$abundance_colnames
                                        test_abundance <- as.matrix(test_abundance)
                                        rownames(test_abundance) <- paste0("test_", seq_len(nrow(test_abundance)))
                                        
                                        # Handle empty samples
                                        row_sums <- rowSums(test_abundance)
                                        empty_idx <- which(row_sums == 0)
                                        valid_idx <- which(row_sums > 0)
                                        
                                        # Initialize predictions
                                        n_test <- nrow(test_abundance)
                                        predictions <- rep(NA_real_, n_test)
                                        
                                        # Process valid (non-empty) samples
                                        if (length(valid_idx) > 0) {
                                          valid_abundance <- test_abundance[valid_idx, , drop = FALSE]
                                          
                                          # Create covariates for valid samples
                                          valid_covariates <- data.frame(
                                            intercept = rep(1, nrow(valid_abundance)),
                                            row.names = rownames(valid_abundance)
                                          )
                                          
                                          # Prepare test data using same offset scheme
                                          test_pln_data <- prepare_data(
                                            counts = valid_abundance,
                                            covariates = valid_covariates, 
                                            offset = self$model$params$offset_scheme
                                          )
                                          
                                          # Get predictions from both models
                                          glmnet_pred <- private$.predict_glmnet(test_pln_data)
                                          pln_pred <- private$.predict_pln(test_pln_data)
                                          
                                          # Combine predictions
                                          w <- self$model$params$weight_glmnet
                                          
                                          if (self$model$glmnet_success && self$model$pln_success) {
                                            # Both models available - use weighted combination
                                            predictions[valid_idx] <- w * glmnet_pred + (1 - w) * pln_pred
                                            
                                          } else if (self$model$glmnet_success) {
                                            # Only glmnet available
                                            predictions[valid_idx] <- glmnet_pred
                                            
                                          } else if (self$model$pln_success) {
                                            # Only PLN available  
                                            predictions[valid_idx] <- pln_pred
                                            
                                          } else {
                                            # Neither model available - use fallback
                                            predictions[valid_idx] <- mean(test_pln_data$Abundance[, self$model$target_name])
                                          }
                                        }
                                        
                                        # Handle empty samples
                                        if (length(empty_idx) > 0) {
                                          predictions[empty_idx] <- 0
                                        }
                                        
                                        # Fill any remaining NAs with fallback
                                        if (any(is.na(predictions))) {
                                          fallback_value <- mean(predictions, na.rm = TRUE)
                                          if (is.na(fallback_value)) fallback_value <- 0
                                          predictions[is.na(predictions)] <- fallback_value
                                        }
                                        
                                        return(list(response = predictions))
                                      },
                                      
                                      #' @description Get cv.glmnet predictions
                                      .predict_glmnet = function(test_pln_data) {
                                        if (!self$model$glmnet_success) {
                                          return(rep(0, nrow(test_pln_data$Abundance)))
                                        }
                                        
                                        tryCatch({
                                          X_test <- test_pln_data$Abundance[, self$model$feature_names, drop = FALSE]
                                          
                                          # Use selected features if specified
                                          if (self$model$params$use_selected_features && length(self$model$selected_features) > 0) {
                                            # Focus on selected features for prediction
                                            X_test <- X_test[, self$model$selected_features, drop = FALSE]
                                            
                                            # Refit a simpler model with selected features only
                                            y_train <- self$model$pln_data$Abundance[, self$model$target_name]
                                            X_train <- self$model$pln_data$Abundance[, self$model$selected_features, drop = FALSE]
                                            
                                            simple_model <- cv.glmnet(
                                              x = X_train,
                                              y = y_train,
                                              family = self$model$params$family_glmnet,
                                              alpha = self$model$params$alpha,
                                              nfolds = self$model$params$nfolds_glmnet
                                            )
                                            
                                            pred <- predict(
                                              simple_model,
                                              newx = X_test,
                                              s = simple_model[[self$model$params$lambda_choice]],
                                              type = "response"
                                            )[, 1]
                                          } else {
                                            # Use full model
                                            pred <- predict(
                                              self$model$glmnet_model,
                                              newx = X_test,
                                              s = self$model$glmnet_model[[self$model$params$lambda_choice]],
                                              type = "response"
                                            )[, 1]
                                          }
                                          
                                          return(pred)
                                          
                                        }, error = function(e) {
                                          # Fallback to mean prediction
                                          return(rep(mean(self$model$pln_data$Abundance[, self$model$target_name]), 
                                                     nrow(test_pln_data$Abundance)))
                                        })
                                      },
                                      
                                      #' @description Get PLN predictions using conditional prediction
                                      .predict_pln = function(test_pln_data) {
                                        if (!self$model$pln_success) {
                                          return(rep(0, nrow(test_pln_data$Abundance)))
                                        }
                                        
                                        tryCatch({
                                          # Create newdata for PLN prediction (covariates)
                                          newdata <- data.frame(
                                            "(Intercept)" = rep(1, nrow(test_pln_data$Abundance)),
                                            row.names = rownames(test_pln_data$Abundance),
                                            check.names = FALSE
                                          )
                                          
                                          # Conditioning data: observed feature species
                                          # Use selected features if available and specified
                                          if (self$model$params$use_selected_features && length(self$model$selected_features) > 0) {
                                            cond_data <- test_pln_data$Abundance[, self$model$selected_features, drop = FALSE]
                                          } else {
                                            cond_data <- test_pln_data$Abundance[, self$model$feature_names, drop = FALSE]
                                          }
                                          
                                          # Make conditional prediction using real predict_cond method
                                          pred_matrix <- predict_cond(
                                            self$model$pln_model,
                                            newdata = newdata,
                                            cond_responses = cond_data,
                                            type = "response"
                                          )
                                          
                                          # Extract target species prediction (should be first column)
                                          return(pred_matrix[, 1])
                                          
                                        }, error = function(e) {
                                          # Fallback prediction using fitted values
                                          return(rep(mean(fitted(self$model$pln_model)[, self$model$target_name]), 
                                                     nrow(test_pln_data$Abundance)))
                                        })
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
  #glmnet_learner,
  #mlr3::LearnerRegrFeatureless$new(),
  #LearnerRegrPLN$new(),
  LearnerRegrPLNGlmnet$new()
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

#learner_id mean_deviance sd_deviance n_iterations
#<char>         <num>       <num>        <int>
#  1:   regr.cv_glmnet     0.7393526  0.03475649            5
#2: regr.featureless     2.3341718  0.10931648            5
#3:         regr.pln     0.7019022  0.04001535            5
