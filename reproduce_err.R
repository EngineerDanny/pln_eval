library(PLNmodels)
library(data.table)
library(glmnet)
set.seed(1)
data("barents", package = "PLNmodels")
task.dt <- as.data.table(t(barents$Abundance))
task.dt[, names(task.dt) := lapply(.SD, log1p)]
Y <- as.matrix(task.dt)
sample_names <- paste0("s", seq_len(nrow(Y)))
rownames(Y) <- sample_names
prepared_data <- prepare_data(counts = Y, 
                              covariates = data.frame(row.names = sample_names))
input_cols <- colnames(Y)[1:(ncol(Y)-1)]
output_col <- colnames(Y)[ncol(Y)]
fit <- PLNnetwork(
  Abundance ~ 1,
  data    = prepared_data,
  penalties = 0.1, 
  control = PLNnetwork_param(
    backend = "torch"
  ) 
)