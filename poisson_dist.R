
library(PLNmodels)   # ≥ 1.1.0
library(glmnet)      # ≥ 4.1‑x

# ===================================================================
# 1.  LASSO‑INITIALISED PLN  (single response)
# ===================================================================
# Fits one Poisson‑ or Gaussian‑Lasso on C ~ A + B to obtain a sparse
# coefficient vector b_init (length = 2).  Then hands that vector to
# PLN as the starting value for B.
# -------------------------------------------------------------------

lasso_init_single <- function(Y, X, offset = NULL,
                              family = c("poisson", "gaussian"),
                              lambda_choice = "lambda.1se",
                              intercept = TRUE) {
  family <- match.arg(family)
  
  ## run cross‑validated glmnet to pick a sensible λ
  cvfit <- cv.glmnet(x = as.matrix(X),
                     y = if (family == "gaussian") log1p(Y) else Y,
                     family = family,
                     alpha  = 1,              # pure Lasso
                     offset = offset,
                     intercept = intercept)
  
  ## extract coefficients at requested λ
  s <- if (is.character(lambda_choice)) cvfit[[lambda_choice]] else lambda_choice
  beta_hat <- as.numeric(coef(cvfit, s = s))[-1L]   # drop intercept
  beta_hat
}

# --------------------------
# Example minimal workflow  |
# --------------------------
# Suppose `counts` is an n × 3 data.frame named A, B, C
counts  <- data.frame(A = rpois(30, 10),
                       B = rpois(30, 15), 
                      C = rpois(30, 20))
offsetC <- log(rowSums(counts))  # (optional) library size offset
#
X_mat   <- as.matrix(counts[, c("A", "B")])
Y_vec   <- counts$C
#
# ## 1‑a  Get sparse start
b0 <- lasso_init_single(Y_vec, X_mat, offset = offsetC,
                         family = "poisson", lambda_choice = "lambda.1se",
                         intercept = FALSE)
#
# ## 1‑b  Wrap into PLN
ctrl <- PLN_param()
ctrl$inception <- list(B = matrix(b0, nrow = ncol(X_mat)))
#
fit_LassoInit <- PLN(C ~ A + B + offset(offsetC),
                      data    = counts,
                      control = ctrl)

fit_LassoInit
beta_final <- coef(fit_LassoInit)  
beta_final
