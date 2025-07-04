library(PLNmodels)
library(corrplot)
library(glassoFast)
library(igraph)
library(data.table)

set.seed(1)


# Mean Poisson deviance function
meanpoissondeviance <- function(y_true, y_pred) {
  deviance_terms <- ifelse(y_true == 0, 
                           2 * y_pred,
                           2 * (y_true * log(y_true / y_pred) - (y_true - y_pred)))
  return(mean(deviance_terms))
}

# Dataset processing from second code
task.dt <- data.table::fread("~/Projects/pln_eval/data/HMPv13_filtered.csv")
task.dt <- task.dt[1:50, 1:100]
taxa_columns <- setdiff(names(task.dt), "Group_ID")
task.dt[, (taxa_columns) := lapply(.SD, function(x) log1p(x)), .SDcols = taxa_columns]
new_column_names <- paste0("Taxa", taxa_columns)
setnames(task.dt, old = taxa_columns, new = new_column_names)

Y <- as.matrix(task.dt[, .SD, .SDcols = new_column_names])
sample_names <- paste0("s", seq_len(nrow(Y)))
rownames(Y) <- sample_names

prepared_data <- prepare_data(counts = Y, 
                              covariates = data.frame(row.names = sample_names))

# PLNPCA setup similar to first code
nb_cores <- 4
options(future.fork.enable = TRUE)
future::plan("multicore", workers = nb_cores)

n <- nrow(Y) 
p <- ncol(Y)
R_max <- min(n - 1, floor(0.3 * p))

# Fit PLNPCA models
myPLNPCAs <- PLNPCA(Abundance ~ 1, data = prepared_data, ranks = 1:R_max)
plot(myPLNPCAs)
myPLNPCA <- getBestModel(myPLNPCAs)

# Extract covariance matrix and apply glassoFast
S <- sigma(myPLNPCA)
out <- glassoFast(S, rho = 0.1)
Omega_hat <- out$wi

# Build a "fixed-covariance" PLN object that carries Ω̂
ctrl <- PLN_param(covariance = "fixed", Omega = Omega_hat)
pFix <- PLN(Abundance ~ 1, data = prepared_data, control = ctrl)

# Prediction setup following second code approach
input_cols <- colnames(Y)[1:(ncol(Y)-1)]
output_col <- colnames(Y)[ncol(Y)]
new_idx <- tail(sample_names, 5)
Yc <- Y[new_idx, input_cols]
newdat <- data.frame(row.names = new_idx)

# Conditional prediction
pred_output <- predict_cond(pFix, 
                            newdata = newdat, 
                            cond_responses = Yc,
                            type = "response")
print(pred_output)

 # Error calculation from second code
actual_values <- Y[new_idx, output_col]
predicted_values <- as.vector(pred_output)

# Calculate the loss metrics
mpd <- meanpoissondeviance(actual_values, predicted_values)
mae <- mean(abs(actual_values - predicted_values))
rmse <- sqrt(mean((actual_values - predicted_values)^2))

print(paste("Mean Poisson Deviance:", mpd))
print(paste("MAE:", mae))
print(paste("RMSE:", rmse))
Omega_hat
