library(PLNmodels)
library(data.table)
library(glmnet)
set.seed(1)

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



# Step 1: Get principled lambda from cv.glmnet
# Use the same prediction setup as your PLNnetwork
input_cols <- colnames(Y)[1:(ncol(Y)-1)]
output_col <- colnames(Y)[ncol(Y)]

cv_fit <- cv.glmnet(Y[, input_cols], Y[, output_col], 
                    family = "poisson", alpha = 1)

# Step 2: Create focused grid around cv-selected lambda
center_lambda <- cv_fit$lambda.min
print(paste("CV-selected lambda:", round(center_lambda, 4)))
custom_penalties <- center_lambda * c(5, 10, 100)


# 1-step dense PLN fitted with the torch backend
inception <- PLN(Abundance ~ 1,
                 data    = prepared_data,
                 control = PLN_param(backend = "torch", trace = 0))

fit <- PLNnetwork(
  Abundance ~ 1,
  data    = prepared_data,
  penalties = 0.01, 
  #penalties = custom_penalties,
  control = PLNnetwork_param(
#    backend = "torch",
    inception = inception
    ) 
)

fit <- PLNnetwork(Abundance ~ 1, data = prepared_data,
                  #penalties = custom_penalties,
                  #penalties = 10, 
                  control = PLNnetwork_param(backend = "nlopt"))
#fit <- PLN(Abundance ~ 1, data = prepared_data)
net_model <- fit

#net_model <- getModel(fit, lambda0)
net_model <- getBestModel(fit)

new_idx <- tail(sample_names, 5)
input_cols <- colnames(Y)[1:(ncol(Y)-1)]
output_col <- colnames(Y)[ncol(Y)]
Yc <- Y[new_idx, input_cols]
newdat <- data.frame(row.names = new_idx)


full_sigma <- solve(net_model$model_par$Omega)  # Full covariance
input_cols <- colnames(Y)[1:(ncol(Y)-1)]
conditioning_sigma <- full_sigma[input_cols, input_cols]
condition_number <- kappa(conditioning_sigma)

print(paste("Conditioning submatrix condition number:", condition_number))

if(condition_number > 1e12) {
  print("Submatrix is nearly singular!")
}

pred_output <- predict_cond(net_model,
                            newdata = newdat,
                            cond_responses = Yc,
                            type = "response")
print(pred_output)



# Get the actual values for comparison
actual_values <- Y[new_idx, output_col]

# Extract predictions (assuming pred_output is a matrix/vector)
predicted_values <- as.vector(pred_output)

# Calculate mean Poisson deviance
meanpoissondeviance <- function(y_true, y_pred) {
  deviance_terms <- ifelse(y_true == 0, 
                           2 * y_pred,
                           2 * (y_true * log(y_true / y_pred) - (y_true - y_pred)))
  return(mean(deviance_terms))
}

# Calculate the loss
mpd <- meanpoissondeviance(actual_values, predicted_values)
mae <- mean(abs(actual_values - predicted_values))
rmse <- sqrt(mean((actual_values - predicted_values)^2))

print(paste("Mean Poisson Deviance:", mpd))
# "Mean Poisson Deviance: 0.857960467950136"
print(paste("MAE:", mae))
print(paste("RMSE:", rmse))