source("/projects/genomic-ml/da2343/PLN/pln_eval/load_source.R")
set.seed(42)
n <- 30
p <- 5
counts <- matrix(rpois(n * p, lambda = 10), nrow = n, ncol = p)
rownames(counts) <- paste0("S", seq_len(n))
colnames(counts) <- paste0("Taxa", seq_len(p))
task <- mlr3::TaskRegr$new("test", as.data.frame(counts), target = "Taxa1")

pln <- LearnerRegrPLN$new()
pln$train(task)
pred <- pln$predict(task)
stopifnot(length(pred$response) == n, !anyNA(pred$response))
cat("LearnerRegrPLN: OK\n")

net <- LearnerRegrPLNnetwork$new()
net$train(task)
pred <- net$predict(task)
stopifnot(length(pred$response) == n, !anyNA(pred$response))
cat("LearnerRegrPLNnetwork: OK\n")

pca <- LearnerRegrPLNPCA$new()
pca$train(task)
pred <- pca$predict(task)
stopifnot(length(pred$response) == n, !anyNA(pred$response))
cat("LearnerRegrPLNPCA: OK\n")

score <- pln$predict(task)$score(poisson_measure)
stopifnot(is.finite(score), score >= 0)
cat("MeasurePoissonDeviance: OK\n")
cat("All tests passed.\n")
