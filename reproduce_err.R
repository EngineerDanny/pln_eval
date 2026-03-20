library(mlr3)
library(PLNmodels)
source("/projects/genomic-ml/da2343/PLN/pln_eval/load_source.R")

future::plan("sequential")

raw <- read.table(
  gzfile("/projects/genomic-ml/da2343/PLN/pln_eval/data/family/microbiomehd/crc_baxter_family.tsv.gz"),
  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
)
counts <- log1p(t(as.matrix(raw)))
counts <- counts[rowSums(counts) > 0, colSums(counts) > 0]
dt <- as.data.frame(counts)
colnames(dt) <- paste0("Taxa", colnames(dt))

target <- colnames(dt)[1]
cat("Target:", target, "\n")
task <- mlr3::TaskRegr$new("test", dt, target = target)

# Run exactly as batchtools would: sequential resample, 3-fold CV
learner <- LearnerRegrPLN$new()
cat("Backend:", learner$param_set$values$backend, "\n")

# No encapsulate — let the raw error surface
cv <- mlr3::rsmp("cv", folds = 3)
cv$instantiate(task)
fold_ids <- cv$train_set(1)

train_task <- task$clone()
train_task$filter(fold_ids)
cat("Training on", train_task$nrow, "samples,", length(train_task$feature_names), "features\n\n")

# Test resample WITHOUT encapsulate
learner <- LearnerRegrPLN$new()
rr <- mlr3::resample(task, learner, cv)
scores <- rr$score(poisson_measure)
cat("\nPer-fold deviance (no encapsulate):\n")
print(scores[, .(iteration, regr.poisson_deviance)])
