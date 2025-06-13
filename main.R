library(PLNmodels)
library(caret)
library(ggplot2)
library(reshape2)
library(MASS)
library(glmnet)

compare_models <- function(abundance_data, k = 5, seed = 42) {
  abundance_matrix <- abundance_data
  n_samples <- nrow(abundance_matrix)
  n_species <- ncol(abundance_matrix)
  target_species <- 1
  
  cat("=== Cross-Validation Setup ===\n")
  cat("Data:", n_samples, "samples x", n_species, "species\n")
  cat("Target species:", colnames(abundance_matrix)[target_species], "(column", target_species, ")\n")
  cat("Target range:", range(abundance_matrix[, target_species]), ", zeros:", round(mean(abundance_matrix[, target_species] == 0) * 100, 1), "%\n")
  cat("Using", length(2:n_species), "predictor species\n\n")
  
  predictor_species <- 2:n_species
  
  k_folds <- createFolds(1:n_samples, k = k, list = TRUE)
  
  pln_deviance <- vector("list", length(k_folds))
  pln_cond_deviance <- vector("list", length(k_folds))
  glmnet_deviance <- vector("list", length(k_folds))
  baseline_deviance <- vector("list", length(k_folds))
  
  #set.seed(seed)
  
  for(i in seq_along(k_folds)) {
    cat("Processing fold", i, "...\n")
    
    test_index <- k_folds[[i]]
    train_index <- setdiff(1:n_samples, test_index)
    
    train_abundance <- abundance_matrix[train_index, ]
    test_abundance <- abundance_matrix[test_index, ]
    
    test_abundance_predictors <- test_abundance[, predictor_species, drop = FALSE]
    
    target_species_mean <- mean(train_abundance[, 1])
    n_test <- nrow(test_abundance)
    baseline_predictions <- rep(target_species_mean, n_test)
    baseline_predictions <- pmax(baseline_predictions, 1e-10)
    
    poisson_deviance_species <- function(y_obs, y_pred) {
      y_log_term <- ifelse(y_obs == 0, 0, y_obs * log(y_obs / pmax(y_pred, 1e-10)))
      2 * (y_log_term - (y_obs - y_pred))
    }
    
    true_target_species <- test_abundance[, 1]
    
    baseline_deviance[[i]] <- mean(poisson_deviance_species(true_target_species, baseline_predictions))
    
    # Prepare data with proper row names
    train_covariates <- data.frame(Intercept = rep(1, nrow(train_abundance)))
    rownames(train_covariates) <- rownames(train_abundance)
    train_data <- prepare_data(train_abundance, train_covariates)
    
    test_covariates <- data.frame(Intercept = rep(1, nrow(test_abundance)))
    rownames(test_covariates) <- rownames(test_abundance)
    test_data <- prepare_data(test_abundance, test_covariates)
    
    # PLN model
    myPLN <- PLN(Abundance ~ 1, data = train_data)
    pln_predictions_full <- predict(myPLN, newdata = test_data, type = "response")
    
    n_test_after_prepare <- nrow(test_data$Abundance)
    if(n_test_after_prepare != n_test) {
      cat("    WARNING: Sample count mismatch after prepare_data!\n")
    }
    
    true_target_species_matched <- test_abundance[1:n_test_after_prepare, 1]
    pln_pred_target_species <- pln_predictions_full[, 1]
    pln_pred_target_species <- pmax(pln_pred_target_species, 1e-10)
    pln_pred_target_species <- pmin(pln_pred_target_species, 1e6)
    pln_deviance[[i]] <- mean(poisson_deviance_species(true_target_species_matched, pln_pred_target_species))
    
    # PLN conditional prediction
    newX <- data.frame("(Intercept)" = rep(1, n_test), check.names = FALSE)
    rownames(newX) <- rownames(test_abundance)
    
    test_abundance_predictors_df <- as.data.frame(test_abundance_predictors)
    predictor_cols_in_order <- colnames(abundance_matrix)[predictor_species]
    test_abundance_predictors_df <- test_abundance_predictors_df[, predictor_cols_in_order, drop = FALSE]
    
    pln_cond_predictions <- predict_cond(myPLN, newX, test_abundance_predictors_df, type = "response")
    
    # For conditional prediction, we predict only the target species, so it should be column 1
    target_col_in_cond <- 1
    
    pln_cond_pred_target_species <- pln_cond_predictions[, target_col_in_cond]
    pln_cond_pred_target_species <- pmax(pln_cond_pred_target_species, 1e-10)
    pln_cond_pred_target_species <- pmin(pln_cond_pred_target_species, 1e6)
    pln_cond_deviance[[i]] <- mean(poisson_deviance_species(true_target_species, pln_cond_pred_target_species))
    
    # GLMNET model
    X_train <- as.matrix(train_abundance[, predictor_species, drop = FALSE])
    Y_train_target <- train_abundance[, 1]
    X_test <- as.matrix(test_abundance[, predictor_species, drop = FALSE])
    
    # Add more regularization to prevent overfitting
    glmnet_model <- cv.glmnet(X_train, Y_train_target, 
                              #family = "poisson",
                              family = "gaussian",
                              type.measure = "deviance",
                              alpha = 1)
    
    # Use lambda.1se for more conservative regularization
    lambda_choice <- ifelse(length(glmnet_model$lambda.1se) > 0, "lambda.1se", "lambda.min")
    
    lambda_star <- glmnet_model$lambda.min
    glmnet_predictions <- predict(glmnet_model, newx = X_test, s = lambda_star, type = "response")
    
    cat("Lambda min:", glmnet_model$lambda.min, "Lambda 1se:", glmnet_model$lambda.1se, "\n")
    coefs <- coef(glmnet_model, s = lambda_star)
    cat("Non-zero coefficients:", sum(abs(coefs) > 1e-8), "\n")
    
    glmnet_pred_target_species <- as.vector(glmnet_predictions)
    glmnet_pred_target_species <- pmax(glmnet_pred_target_species, 1e-10)
    glmnet_pred_target_species <- pmin(glmnet_pred_target_species, 1e6)
    glmnet_deviance[[i]] <- mean(poisson_deviance_species(true_target_species, glmnet_pred_target_species))
    
    cat("  Fold", i, "results: Baseline =", round(baseline_deviance[[i]], 3), 
        ", PLN =", round(pln_deviance[[i]], 3), 
        ", PLN_COND =", round(pln_cond_deviance[[i]], 3), 
        ", GLMNET =", round(glmnet_deviance[[i]], 3), "\n")
  }
  
  target_name <- colnames(abundance_matrix)[1]
  
  cat("\n=== Final Results - DEVIANCE (Target:", target_name, ") ===\n")
  cat("BASELINE - Mean Deviance:", round(mean(sapply(baseline_deviance, mean)), 4), 
      "± SD:", round(sd(sapply(baseline_deviance, mean)), 4), "\n")
  cat("PLN - Mean Deviance:", round(mean(sapply(pln_deviance, mean)), 4), 
      "± SD:", round(sd(sapply(pln_deviance, mean)), 4), "\n")
  cat("PLN_COND - Mean Deviance:", round(mean(sapply(pln_cond_deviance, mean)), 4), 
      "± SD:", round(sd(sapply(pln_cond_deviance, mean)), 4), "\n")
  cat("GLMNET - Mean Deviance:", round(mean(sapply(glmnet_deviance, mean)), 4), 
      "± SD:", round(sd(sapply(glmnet_deviance, mean)), 4), "\n")
  
  fold_dev_baseline <- sapply(baseline_deviance, mean)
  fold_dev_pln <- sapply(pln_deviance, mean)
  fold_dev_pln_cond <- sapply(pln_cond_deviance, mean)
  fold_dev_glmnet <- sapply(glmnet_deviance, mean)
  
  cat("\n=== Model Comparisons vs Baseline - DEVIANCE (p-values) ===\n")
  cat("BASELINE vs PLN:", round(t.test(fold_dev_baseline, fold_dev_pln, paired = TRUE)$p.value, 4), "\n")
  cat("BASELINE vs PLN_COND:", round(t.test(fold_dev_baseline, fold_dev_pln_cond, paired = TRUE)$p.value, 4), "\n")
  cat("BASELINE vs GLMNET:", round(t.test(fold_dev_baseline, fold_dev_glmnet, paired = TRUE)$p.value, 4), "\n")
  
  cat("\n=== Model Comparisons - DEVIANCE (p-values) ===\n")
  cat("PLN vs PLN_COND:", round(t.test(fold_dev_pln, fold_dev_pln_cond, paired = TRUE)$p.value, 4), "\n")
  cat("PLN vs GLMNET:", round(t.test(fold_dev_pln, fold_dev_glmnet, paired = TRUE)$p.value, 4), "\n")
  cat("PLN_COND vs GLMNET:", round(t.test(fold_dev_pln_cond, fold_dev_glmnet, paired = TRUE)$p.value, 4), "\n")
  
  cat("\n=== Fold-by-fold DEVIANCE ===\n")
  results_df_dev <- data.frame(
    Fold = 1:length(k_folds),
    BASELINE = round(fold_dev_baseline, 4),
    PLN = round(fold_dev_pln, 4),
    PLN_COND = round(fold_dev_pln_cond, 4),
    GLMNET = round(fold_dev_glmnet, 4)
  )
  print(results_df_dev)
  
  cat("\n=== Improvement over Baseline - DEVIANCE ===\n")
  cat("PLN improvement:", round((mean(fold_dev_baseline) - mean(fold_dev_pln)) / mean(fold_dev_baseline) * 100, 2), "%\n")
  cat("PLN_COND improvement:", round((mean(fold_dev_baseline) - mean(fold_dev_pln_cond)) / mean(fold_dev_baseline) * 100, 2), "%\n")
  cat("GLMNET improvement:", round((mean(fold_dev_baseline) - mean(fold_dev_glmnet)) / mean(fold_dev_baseline) * 100, 2), "%\n")
}

filter_sparse_species <- function(abundance_data, min_mean = 0.5, min_prevalence = 0.5, transform = FALSE, transform_method = "log1p") {
  species_means <- colMeans(abundance_data)
  species_prevalence <- colMeans(abundance_data > 0)
  
  keep_species <- (species_means >= min_mean) & (species_prevalence >= min_prevalence)
  
  cat("Removing", sum(!keep_species), "sparse species:\n")
  cat(names(abundance_data)[!keep_species], "\n")
  cat("Keeping", sum(keep_species), "species for analysis\n")
  
  filtered_data <- abundance_data[, keep_species]
  
  if(transform) {
    cat("Applying", transform_method, "transformation:\n")
    for(sp in colnames(filtered_data)) {
      original_range <- range(filtered_data[, sp])
      filtered_data[, sp] <- log1p(filtered_data[, sp])
      new_range <- range(filtered_data[, sp])
      cat("  ", sp, ": range [", original_range[1], ",", original_range[2], 
          "] -> [", round(new_range[1],2), ",", round(new_range[2],2), "]\n")
    }
  }
  
  return(filtered_data)
}

data(trichoptera)
filtered_data <- filter_sparse_species(trichoptera$Abundance)

generate_pln_data <- function(n_samples = 100, n_species = 500, sparsity = 0.9, seed = 123) {
  #set.seed(seed)
  
  # Higher sparsity = lower mu values = more zeros
  mu_min <- 0.5 - sparsity * 1.5  # ranges from 0.5 (sparsity=0) to -1.0 (sparsity=1)
  mu_max <- 2.5 - sparsity * 1.0  # ranges from 2.5 (sparsity=0) to 1.5 (sparsity=1)
  
  mu <- runif(n_species, mu_min, mu_max)
  
  Sigma <- matrix(0.5, n_species, n_species)
  diag(Sigma) <- 1.0
  
  Z <- MASS::mvrnorm(n_samples, mu = mu, Sigma = Sigma)
  
  abundance_matrix <- matrix(0, n_samples, n_species)
  for(i in 1:n_samples) {
    for(j in 1:n_species) {
      abundance_matrix[i, j] <- rpois(1, exp(Z[i, j]))
    }
  }
  
  colnames(abundance_matrix) <- paste0("Sp", 1:n_species)
  rownames(abundance_matrix) <- paste0("Sample_", 1:n_samples)
  
  as.data.frame(abundance_matrix)
}

synth_data <- generate_pln_data()
# print(synth_data)
# Calculate mean and variance for each species

means <- colMeans(synth_data)
variances <- apply(synth_data, 2, var)
dispersion <- variances / means

print(data.frame(
  Species = names(means),
  Mean = round(means, 2),
  Variance = round(variances, 2), 
  Phi = round(dispersion, 2)
))

filter_by_dispersion <- function(abundance_data, phi_min = 0.0, phi_max = 2) {
  means <- colMeans(abundance_data)
  variances <- apply(abundance_data, 2, var)
  phi <- variances / means
  
  keep_species <- (phi >= phi_min) & (phi <= phi_max)
  
  cat("Species dispersion (φ = variance/mean):\n")
  for(i in 1:length(phi)) {
    status <- ifelse(keep_species[i], "KEEP", "REMOVE")
    cat("  ", names(phi)[i], ": φ =", round(phi[i], 2), "(", status, ")\n")
  }
  
  filtered_data <- abundance_data[, keep_species]
  cat("Kept", sum(keep_species), "of", length(keep_species), "species\n")
  
  return(filtered_data)
}

filtered_abund <- filter_by_dispersion(synth_data)
filtered_abund

compare_models(
  abundance_data = filtered_abund,
  #abundance_data = synth_data,
  k = 5
)

