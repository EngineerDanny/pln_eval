library(data.table)
library(glmnet)
library(ggplot2)
library(MASS)

# Load and prepare data (same as original)
task.dt <- data.table::fread("/projects/genomic-ml/da2343/PLN/pln_eval/data/HMPv13_filtered.csv")
task.dt <- task.dt[1:500]  # Use first 500 samples

taxa_columns <- setdiff(names(task.dt), "Group_ID")
new_column_names <- paste0("Taxa", taxa_columns)
setnames(task.dt, old = taxa_columns, new = new_column_names)

# Poisson deviance function (matches your MLR3 version exactly)
poisson_deviance <- function(y_true, y_pred) {
  eps <- 1e-10
  y_pred <- pmax(y_pred, eps)
  log_term <- ifelse(y_true == 0, 0, y_true * log(y_true / y_pred))
  2 * mean(log_term - (y_true - y_pred))
}

# Set up target analysis
target_taxa <- "Taxa4365684"
cat("=== Analyzing", target_taxa, "===\n")

# Extract target and features
y <- task.dt[[target_taxa]]
feature_names <- new_column_names[new_column_names != target_taxa]
X <- as.matrix(task.dt[, ..feature_names])

# Check target characteristics
cat("Target stats: mean =", round(mean(y), 2), 
    ", zeros =", round(mean(y == 0) * 100, 1), "%\n")
cat("Variance:", round(var(y), 2), ", Dispersion (φ):", round(var(y)/mean(y), 1), "\n")

# Calculate optimal theta for negative binomial
set.seed(42)
optimal_theta <- theta.ml(y, mu = mean(y))
cat("Optimal theta:", round(optimal_theta, 4), "\n")

# Set up 5-fold cross-validation (matching your original setup)
n_folds <- 5
n_samples <- length(y)
fold_indices <- sample(rep(1:n_folds, length.out = n_samples))

# Storage for results
glmnet_deviances <- numeric(n_folds)
featureless_deviances <- numeric(n_folds)
glmnet_success <- logical(n_folds)

cat("\n=== Running 5-fold Cross-Validation ===\n")

# Cross-validation loop
for(fold in 1:n_folds) {
  cat("Processing fold", fold, "/", n_folds, "...")
  
  # Split data (same logic as MLR3 would use)
  test_idx <- which(fold_indices == fold)
  train_idx <- which(fold_indices != fold)
  
  X_train <- X[train_idx, , drop = FALSE]
  y_train <- y[train_idx]
  X_test <- X[test_idx, , drop = FALSE]
  y_test <- y[test_idx]
  
  # 1. Featureless baseline (mean prediction - same as MLR3 featureless)
  mean_pred <- mean(y_train)
  featureless_pred <- rep(mean_pred, length(y_test))
  featureless_deviances[fold] <- poisson_deviance(y_test, featureless_pred)
  
  # 2. GLMNet with negative binomial
  tryCatch({
    # Use cv.glmnet with negative binomial (this works with direct glmnet)
    cv_nb_model <- cv.glmnet(
      x = X_train, 
      y = y_train,
      family = negative.binomial(theta = optimal_theta),
      alpha = 1,  # Lasso (same as your original)
      nfolds = 5,  # Internal CV for lambda selection
      standardize = TRUE,  # Default standardization
      maxit = 1000000  # Prevent convergence issues
    )
    
    # Predict on test set using optimal lambda
    glmnet_pred <- predict(cv_nb_model, newx = X_test, s = "lambda.min", type = "response")
    glmnet_pred <- as.numeric(glmnet_pred)
    
    # Calculate Poisson deviance
    glmnet_deviances[fold] <- poisson_deviance(y_test, glmnet_pred)
    glmnet_success[fold] <- TRUE
    
    # Count selected features
    coefs <- coef(cv_nb_model, s = "lambda.min")
    n_features <- sum(abs(coefs[-1]) > 1e-8)  # Exclude intercept
    cat(" Features:", n_features, "✓\n")
    
  }, error = function(e) {
    cat(" ERROR:", e$message, "\n")
    glmnet_deviances[fold] <- NA
    glmnet_success[fold] <- FALSE
  })
}

# Calculate aggregate results (matching your original format)
successful_folds <- which(glmnet_success)
n_successful <- length(successful_folds)

if(n_successful >= 3) {
  # Calculate means and SDs (using log as in your original)
  mean_glmnet_log <- mean(log(glmnet_deviances[successful_folds]), na.rm = TRUE)
  sd_glmnet_log <- sd(log(glmnet_deviances[successful_folds]), na.rm = TRUE)
  mean_featureless_log <- mean(log(featureless_deviances[successful_folds]), na.rm = TRUE)
  sd_featureless_log <- sd(log(featureless_deviances[successful_folds]), na.rm = TRUE)
  
  # Create results table (matching your original format)
  aggregate_results <- data.table(
    learner_id = c("regr.cv_glmnet", "regr.featureless"),
    train.subsets = c("all", "all"),
    mean_deviance = c(mean_glmnet_log, mean_featureless_log),
    sd_deviance = c(sd_glmnet_log, sd_featureless_log),
    n_iterations = c(n_successful, n_successful)
  )
  
  cat("\n=== Results (log scale, matching original) ===\n")
  print(aggregate_results)
  
  # Statistical test (same as your original)
  t_test_result <- t.test(
    glmnet_deviances[successful_folds],
    featureless_deviances[successful_folds],
    paired = TRUE
  )
  
  cat("\n=== Statistical Test ===\n")
  print(t_test_result)
  
  # Additional summary (raw deviance scale)
  cat("\n=== Raw Deviance Summary ===\n")
  cat("GLMNet NB - Mean:", round(mean(glmnet_deviances[successful_folds]), 4), 
      "± SD:", round(sd(glmnet_deviances[successful_folds]), 4), "\n")
  cat("Featureless - Mean:", round(mean(featureless_deviances[successful_folds]), 4), 
      "± SD:", round(sd(featureless_deviances[successful_folds]), 4), "\n")
  
  improvement <- (mean(featureless_deviances[successful_folds]) - mean(glmnet_deviances[successful_folds])) / 
    mean(featureless_deviances[successful_folds]) * 100
  cat("Improvement:", round(improvement, 2), "%\n")
  cat("Successful folds:", n_successful, "/", n_folds, "\n")
  
} else {
  cat("ERROR: Too few successful folds (", n_successful, "/", n_folds, ")\n")
}

# Target characteristics (matching your original output)
cat("\n=== Target Variable Summary ===\n")
summary(task.dt[[target_taxa]])
cat("Zeros:", mean(task.dt[[target_taxa]] == 0) * 100, "%\n")
cat("Mean:", mean(task.dt[[target_taxa]]), "\n")
cat("Variance:", var(task.dt[[target_taxa]]), "\n")