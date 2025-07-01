compute_PLN_starting_point <- function(Y, X, O, w, s = 0.1, use_glmnet = F) {
  n <- nrow(Y); p <- ncol(Y); d <- ncol(X)
  
  if (use_glmnet && requireNamespace("glmnet", quietly = TRUE)) {
    Y_transformed <- log((1 + Y) / exp(O))
    X_weighted <- sqrt(w) * X
    
    # Initialize coefficient matrix - ensure correct dimensions
    B <- matrix(0, nrow = d, ncol = p)
    
    for (j in 1:p) {
      y_j_weighted <- sqrt(w) * Y_transformed[, j]
      
      cv_fit <- cv.glmnet(
        x = X_weighted,
        y = y_j_weighted,
        family = "gaussian",
        alpha = 1,
        standardize = TRUE,
        intercept = FALSE  # Match original approach
      )
      
      # Extract coefficients correctly - exclude intercept
      coefs <- as.vector(coef(cv_fit, s = "lambda.min"))[-1]  # Remove intercept
      B[, j] <- coefs
    }
    
    # Use weighted prediction for consistency
    fits_pred <- X_weighted %*% B
    residuals_matrix <- Y_transformed - fits_pred / sqrt(w)  # Adjust for weighting
    
  } else {
    fits <- lm.fit(w * X, w * log((1 + Y)/exp(O)))
    B <- matrix(fits$coefficients, d, p)
    residuals_matrix <- matrix(fits$residuals, n, p)
  }
  
  list(
    B = B,
    M = residuals_matrix,
    S = matrix(s, n, p)
  )
}


set.seed(123)
n <- 50; p <- 3; d <- 2

# Create true relationships
B_true <- matrix(c(0.5, -0.3, 0.2, 0.1, -0.4, 0.6), nrow = d, ncol = p)
X <- matrix(rnorm(n * d), nrow = n, ncol = d)
linear_pred <- X %*% B_true
Y <- matrix(rpois(n * p, lambda = exp(linear_pred)), nrow = n, ncol = p)
O <- matrix(0, nrow = n, ncol = p)
w <- rep(1, n)

result <- PLNmodels::compute_PLN_starting_point(Y, X, O, w, s = 0.1)
print("True B:")
print(B_true)
print("Estimated B:")
print(result$B)
