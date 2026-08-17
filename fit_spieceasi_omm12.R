#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(SpiecEasi)
})

base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
proc_dir <- file.path(base_dir, "data", "interaction_ground_truth", "omm12", "processed")
input_path <- file.path(proc_dir, "community_absabundance_in_vitro.tsv.gz")
output_path <- file.path(proc_dir, "benchmark_outputs", "SPIECEASI_predictions.tsv.gz")

abundance <- fread(cmd = paste("gzip -cd", shQuote(input_path)))
X <- as.matrix(abundance[, -1, with = FALSE])
rownames(X) <- abundance[[1]]

set.seed(426001)
fit <- spiec.easi(
  X,
  method = "glasso",
  nlambda = 50,
  lambda.min.ratio = 0.05,
  verbose = FALSE
)

opt_index <- getOptInd(fit)
lambda_path <- fit$lambda
cat(
  "Selected lambda:", getOptLambda(fit),
  "| minimum lambda:", min(lambda_path),
  "| selected minimum:", opt_index == length(lambda_path), "\n"
)

adjacency <- as.matrix(getRefit(fit))
precision <- as.matrix(getOptiCov(fit))
partial_cor <- -precision / sqrt(outer(diag(precision), diag(precision)))
diag(partial_cor) <- 0

edge_idx <- which(upper.tri(adjacency) & adjacency != 0, arr.ind = TRUE)
edges <- data.table(
  taxon_1 = colnames(X)[edge_idx[, 1]],
  taxon_2 = colnames(X)[edge_idx[, 2]],
  sign = ifelse(partial_cor[edge_idx] > 0, 1L, -1L),
  weight = partial_cor[edge_idx]
)
fwrite(edges, output_path, sep = "\t", compress = "gzip")
cat("Saved", nrow(edges), "SPIEC-EASI edges to", output_path, "\n")
