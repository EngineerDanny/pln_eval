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
library(PLNmodels)
library(corrplot)
library(glassoFast)
library(igraph)


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
                                  pln_model <- PLN(Abundance ~ 1, 
                                                   data = pln_data, 
                                                   control = control_params)
                                  
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

LearnerRegrCVPLN <- R6::R6Class("LearnerRegrCVPLN",
                                inherit = LearnerRegr,
                                public = list(
                                  initialize = function() {
                                    ps = ps(
                                      alpha = p_dbl(0, 1, default = 1, tags = "train"),
                                      lambda_multipliers = p_uty(default = c(0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5), tags = "train"),
                                      max_condition_number = p_dbl(1, 1e15, default = 1e12, tags = "train"),
                                      covariance = p_fct(c("full", "diagonal", "spherical"), default = "full", tags = "train"),
                                      trace = p_int(0, 2, default = 0, tags = "train"),
                                      backend = p_fct(c("nlopt", "torch"), default = "nlopt", tags = "train"),
                                      offset_scheme = p_fct(c("TSS", "CSS", "RLE", "GMPR", "none"), default = "TSS", tags = "train")
                                    )
                                    ps$values = list(
                                      alpha = 1,
                                      lambda_multipliers = c(0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5),
                                      max_condition_number = 1e12,
                                      covariance = "full", 
                                      trace = 0, 
                                      backend = "nlopt", 
                                      offset_scheme = "TSS"
                                    )
                                    super$initialize(
                                      id = "regr.cv_pln",
                                      param_set = ps,
                                      feature_types = c("integer", "numeric"),
                                      label = "Poisson Log-Normal Model with CV Lambda Selection",
                                      packages = c("PLNmodels", "glmnet")
                                    )
                                  }
                                ),
                                private = list(
                                  .train = function(task) {
                                    pv <- self$param_set$get_values(tags = "train")
                                    target_name <- task$target_names
                                    feature_names <- task$feature_names
                                    
                                    # Get data and convert to matrix
                                    task_data <- task$data(cols = c(target_name, feature_names))
                                    abundance_matrix <- data.matrix(task_data, rownames.force = TRUE)
                                    
                                    # Prepare covariates (intercept only)
                                    covariates <- data.frame(row.names = rownames(abundance_matrix))
                                    
                                    # Prepare PLN data
                                    pln_data <- prepare_data(
                                      counts = abundance_matrix, 
                                      covariates = covariates, 
                                      offset = pv$offset_scheme
                                    )
                                    
                                    inception <- PLN(Abundance ~ 1,
                                                     data = pln_data,
                                                     control = PLN_param(backend = "torch", 
                                                                         trace = 0))
                                    
                                  
                                    net_fit <- PLNnetwork(
                                      Abundance ~ 1,
                                      data = pln_data,
                                      #penalties = center_lambda/10,
                                      penalties = 0.1,
                                      control = PLNnetwork_param(
                                        inception = inception
                                      ) 
                                    )
                                    
                                    # Get best model
                                    best_model <- getBestModel(net_fit)
                                    
                                    # Check condition number for numerical stability
                                    full_sigma <- solve(best_model$model_par$Omega)
                                    conditioning_sigma <- full_sigma[input_cols, input_cols]
                                    condition_number <- kappa(conditioning_sigma)
                                    
                                    if (pv$trace > 0) {
                                      cat("Selected penalty:", round(best_model$penalty, 4), "\n")
                                      cat("Conditioning submatrix condition number:", round(condition_number, 2), "\n")
                                    }
                                    
                                    if (condition_number > pv$max_condition_number) {
                                      warning(paste("Conditioning submatrix condition number:", 
                                                    round(condition_number), "exceeds threshold of", 
                                                    pv$max_condition_number, "- predictions may be unstable"))
                                    }
                                    
                                    # Store model and metadata
                                    self$model <- list(
                                      network_model = best_model,
                                      target_name = target_name, 
                                      feature_names = feature_names,
                                      center_lambda = center_lambda,
                                      selected_penalty = best_model$penalty,
                                      condition_number = condition_number,
                                      network_density = best_model$density
                                    )
                                    
                                    invisible(self$model)
                                  },
                                  .predict = function(task) {
                                    target_name <- self$model$target_name
                                    feature_names <- self$model$feature_names
                                    
                                    # Get prediction data
                                    task_data <- task$data(cols = c(target_name, feature_names))
                                    abundance_matrix <- data.matrix(task_data, rownames.force = TRUE)
                                    
                                    # Prepare covariates for prediction (empty data.frame with row names)
                                    covariates <- data.frame(row.names = rownames(abundance_matrix))
                                    
                                    # Get conditioning data (input features)
                                    conditioning_data <- abundance_matrix[, feature_names, drop = FALSE]
                                    
                                    # Make predictions using predict_cond
                                    pred <- predict_cond(
                                      self$model$network_model, 
                                      newdata = covariates, 
                                      cond_responses = conditioning_data, 
                                      type = "response"
                                    )
                                    
                                    # Extract predictions (assuming single target column)
                                    pred_values <- as.vector(pred[, 1])
                                    
                                    list(response = pred_values)
                                  }
                                )
)

LearnerRegrPLNPCA <- R6::R6Class("LearnerRegrPLNPCA",
                                 inherit = LearnerRegr,
                                 public = list(
                                   initialize = function() {
                                     ps = ps(
                                       covariance = p_fct(c("full", "diagonal", "spherical"), default = "full", tags = "train"),
                                       trace = p_int(0, 2, default = 0, tags = "train"),
                                       backend = p_fct(c("nlopt", "torch"), default = "nlopt", tags = "train"),
                                       offset_scheme = p_fct(c("TSS", "CSS", "RLE", "GMPR", "none"), default = "TSS", tags = "train"),
                                       rank_max = p_int(2L, default = 5L, tags = "train"),
                                       criterion = p_fct(c("BIC", "ICL"), default = "BIC", tags = "train"),
                                       rho = p_dbl(0.001, default = 0.1, tags = "train")
                                     )
                                     ps$values = list(covariance = "full", trace = 0, backend = "nlopt", 
                                                      offset_scheme = "TSS", rank_max = 5L, criterion = "BIC", rho = 0.1)
                                     super$initialize(
                                       id = "regr.plnpca",
                                       param_set = ps,
                                       feature_types = c("integer", "numeric"),
                                       label = "Poisson Log-Normal PCA Model",
                                       packages = c("PLNmodels", "mlr3learners", "glassoFast")
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
                                     pln_data <- prepare_data(
                                       counts = abundance_matrix, 
                                       covariates = covariates,
                                       offset = pv$offset_scheme %||% "TSS"
                                     )
                                     
                                     # Calculate ranks range based on data dimensions and parameters
                                     n <- nrow(pln_data$Abundance)
                                     p <- ncol(pln_data$Abundance)
                                     ranks <- 1:min(n - 1, floor(0.3 * p))
                                     
                                     # Fit PLNPCA model family on prepared data
                                     plnpca_family <- PLNPCA(Abundance ~ 1, 
                                                             data = pln_data, 
                                                             ranks = ranks)
                                     
                                     # Extract best model using specified criterion
                                     best_model <- getBestModel(plnpca_family, pv$criterion %||% "BIC")
                                     
                                     # Extract covariance matrix and apply glassoFast
                                     rho_param <- 0.5 
                                     out <- glassoFast::glassoFast(sigma(best_model), rho = rho_param)
                                     Omega_hat <- out$wi
                                     
                                     # Build a "fixed-covariance" PLN object that carries Ω̂
                                     ctrl <- PLN_param(covariance = "fixed", Omega = Omega_hat)
                                     pFix <- PLN(Abundance ~ 1, data = pln_data, control = ctrl)
                                     
                                     # Store model and related info
                                     self$model <- list(
                                       pln_model = pFix,
                                       plnpca_family = plnpca_family,
                                       best_plnpca = best_model,
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
                                       
                                       # Make conditional predictions using predict_cond on fixed PLN model
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
