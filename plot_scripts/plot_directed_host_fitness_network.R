library(data.table)
library(igraph)
library(glmnet)

base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
truth_dir <- file.path(base_dir, "data", "interaction_ground_truth", "host_fitness_2018", "processed")
fig_dir <- file.path(base_dir, "figures", "march26")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv_gz <- function(path) {
  read.delim(gzfile(path), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

prepare_count_like <- function(X, target_depth = 10000L) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  X[!is.finite(X)] <- 0
  X[X < 0] <- 0
  is_integerish <- all(abs(X - round(X)) < 1e-8)
  if (is_integerish) return(round(X))
  rs <- rowSums(X)
  scaled <- matrix(0, nrow(X), ncol(X), dimnames = dimnames(X))
  keep <- rs > 0
  scaled[keep, ] <- round(X[keep, , drop = FALSE] / rs[keep] * target_depth)
  scaled
}

safe_glmnet_coefs <- function(x_train, y_train) {
  p_minus_1 <- ncol(x_train)
  zero_coefs <- rep(0, p_minus_1)
  if (p_minus_1 == 0L || length(unique(y_train)) < 2L) return(zero_coefs)
  x_train <- as.matrix(x_train)
  cv_nfolds <- max(3L, min(5L, nrow(x_train)))
  fit <- tryCatch(cv.glmnet(x_train, y_train, family = "poisson", alpha = 1, nfolds = cv_nfolds), error = function(e) NULL)
  if (is.null(fit)) return(zero_coefs)
  coefs <- tryCatch(as.numeric(coef(fit, s = "lambda.1se")[-1]), error = function(e) zero_coefs)
  if (length(coefs) != p_minus_1) zero_coefs else coefs
}

fit_loto_glmnet_coef_mat <- function(X_train_raw, X_train_log) {
  p <- ncol(X_train_raw)
  coef_mat <- matrix(0, p, p, dimnames = list(colnames(X_train_raw), colnames(X_train_raw)))
  for (j in seq_len(p)) {
    y_train <- X_train_raw[, j]
    x_train <- X_train_log[, -j, drop = FALSE]
    coef_mat[j, -j] <- safe_glmnet_coefs(x_train, y_train)
  }
  coef_mat
}

curvature_for_directed <- function(edge_dt) {
  curv <- rep(0, nrow(edge_dt))
  if (!nrow(edge_dt)) return(curv)
  keys <- paste(pmin(edge_dt$from, edge_dt$to), pmax(edge_dt$from, edge_dt$to), sep = "||")
  for (k in unique(keys)) {
    idx <- which(keys == k)
    if (length(idx) == 2L) {
      ord <- order(paste(edge_dt$from[idx], edge_dt$to[idx], sep = "->"))
      curv[idx[ord[1]]] <- 0.22
      curv[idx[ord[2]]] <- -0.22
    }
  }
  curv
}

plot_frame <- function() {
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], border = "grey45", lwd = 1.6, xpd = NA)
}

taxa <- c("LP", "LB", "AP", "AT", "AO")
coords <- rbind(
  LP = c(-0.9, 0.2),
  AO = c(0.0, 0.95),
  AP = c(1.0, -0.1),
  LB = c(0.05, -0.95),
  AT = c(-1.05, -0.15)
)

truth <- as.data.table(read_tsv_gz(file.path(truth_dir, "truth_directed.tsv.gz")))
truth_edges <- truth[, .(
  from = source_taxon,
  to = target_taxon,
  weight = as.numeric(max_abs_strength),
  edge_sign = ifelse(effect_sign > 0, "positive", "negative")
)]

pln_pred <- as.data.table(read_tsv_gz(file.path(truth_dir, "benchmark_outputs", "PLNnetwork_predictions.tsv.gz")))
pln_edges <- pln_pred[, .(
  from = taxon_1,
  to = taxon_2,
  weight = abs(as.numeric(weight)),
  edge_sign = ifelse(sign > 0, "positive", "negative")
)]

abundance <- as.data.table(read_tsv_gz(file.path(truth_dir, "abundance_matrix.tsv.gz")))
X_counts <- prepare_count_like(as.matrix(abundance[, ..taxa]))
colnames(X_counts) <- taxa
X_log <- log1p(X_counts)
glm_coef <- fit_loto_glmnet_coef_mat(X_counts, X_log)
glm_idx <- which(abs(glm_coef) > 1e-8, arr.ind = TRUE)
glm_idx <- glm_idx[glm_idx[, 1] != glm_idx[, 2], , drop = FALSE]
glm_edges <- data.table(
  from = rownames(glm_coef)[glm_idx[, 1]],
  to = colnames(glm_coef)[glm_idx[, 2]],
  weight = abs(glm_coef[glm_idx]),
  edge_sign = ifelse(glm_coef[glm_idx] > 0, "positive", "negative")
)

plot_directed_graph <- function(edge_dt, main_title, show_legend = FALSE) {
  g <- graph_from_data_frame(edge_dt, directed = TRUE, vertices = data.frame(name = taxa))
  eweights <- E(g)$weight
  edge_widths <- if (length(eweights)) 1 + 5 * (eweights / max(eweights)) else numeric()
  edge_colors <- ifelse(E(g)$edge_sign == "positive", "#1f78b4", "#e31a1c")
  curv <- curvature_for_directed(edge_dt)

  plot(
    g,
    layout = coords[V(g)$name, , drop = FALSE],
    vertex.label = V(g)$name,
    vertex.label.cex = 0.9,
    vertex.label.family = "sans",
    vertex.label.color = "black",
    vertex.size = 24,
    vertex.color = "grey96",
    vertex.frame.color = "grey35",
    edge.width = edge_widths,
    edge.color = edge_colors,
    edge.curved = curv,
    edge.arrow.size = 0.45,
    edge.arrow.width = 1.0,
    main = main_title,
    margin = 0.05
  )
  plot_frame()
  if (show_legend) {
    legend(
      "bottomleft",
      legend = c("Positive", "Negative"),
      col = c("#1f78b4", "#e31a1c"),
      lwd = c(2.5, 2.5),
      bty = "n",
      cex = 0.8,
      inset = c(0.01, 0.02)
    )
  }
}

plot_undirected_graph <- function(edge_dt, main_title) {
  g <- graph_from_data_frame(edge_dt, directed = FALSE, vertices = data.frame(name = taxa))
  eweights <- E(g)$weight
  edge_widths <- if (length(eweights)) 1 + 5 * (eweights / max(eweights)) else numeric()
  edge_colors <- ifelse(E(g)$edge_sign == "positive", "#1f78b4", "#e31a1c")

  plot(
    g,
    layout = coords[V(g)$name, , drop = FALSE],
    vertex.label = V(g)$name,
    vertex.label.cex = 0.9,
    vertex.label.family = "sans",
    vertex.label.color = "black",
    vertex.size = 24,
    vertex.color = "grey96",
    vertex.frame.color = "grey35",
    edge.width = edge_widths,
    edge.color = edge_colors,
    edge.curved = 0.08,
    main = main_title,
    margin = 0.05
  )
  plot_frame()
}

out_png <- file.path(fig_dir, "interaction_network_directed_host_fitness_2018.png")
png(out_png, width = 9.2, height = 3.9, units = "in", res = 300)
op <- par(mfrow = c(1, 3), mar = c(0.3, 0.3, 2.1, 0.3), oma = c(0, 0, 0, 0))
plot_directed_graph(truth_edges, "Truth (directed)", show_legend = TRUE)
plot_undirected_graph(pln_edges, "PLNNetwork (undirected)")
plot_directed_graph(glm_edges, "LOTO glmnet (directed)")
par(op)
dev.off()

cat(out_png, "\n")
