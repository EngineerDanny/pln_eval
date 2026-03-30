library(data.table)
library(igraph)

args <- commandArgs(trailingOnly = TRUE)
base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
truth_dir <- file.path(base_dir, "data", "interaction_ground_truth")
fig_dir <- file.path(base_dir, "figures", "march26")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
date_tag <- format(Sys.Date(), "%Y-%m-%d")

label_map <- c(
  omm12 = "OMM12",
  omm12_keystone_2023 = "OMM12 keystone 2023",
  host_fitness_2018 = "Host fitness 2018",
  butyrate_assembly_2021 = "Butyrate assembly 2021",
  pairinterax = "PairInteraX"
)

dataset <- if (length(args) >= 1) args[1] else "omm12_keystone_2023"
dataset_label <- if (!is.na(label_map[[dataset]])) label_map[[dataset]] else dataset
proc_dir <- file.path(truth_dir, dataset, "processed")

read_tsv_gz <- function(path) {
  read.delim(gzfile(path), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

truth <- as.data.table(read_tsv_gz(file.path(proc_dir, "truth_undirected.tsv.gz")))
pln <- as.data.table(read_tsv_gz(file.path(proc_dir, "benchmark_outputs", "PLNnetwork_predictions.tsv.gz")))
glm <- as.data.table(read_tsv_gz(file.path(proc_dir, "benchmark_outputs", "LOTO_glmnet_CV1se_predictions.tsv.gz")))

truth_edges <- truth[, .(
  from = taxon_1,
  to = taxon_2,
  weight = as.numeric(max_abs_strength),
  edge_sign = ifelse(
    sign_consensus > 0, "positive",
    ifelse(sign_consensus < 0, "negative", "mixed")
  )
)]
pln_edges <- pln[, .(
  from = taxon_1,
  to = taxon_2,
  weight = abs(as.numeric(weight)),
  edge_sign = ifelse(sign > 0, "positive", "negative")
)]
glm_edges <- glm[, .(
  from = taxon_1,
  to = taxon_2,
  weight = abs(as.numeric(weight)),
  edge_sign = ifelse(sign > 0, "positive", "negative")
)]

union_edges <- unique(rbindlist(list(
  truth_edges[, .(from, to)],
  pln_edges[, .(from, to)],
  glm_edges[, .(from, to)]
)))
vertices <- sort(unique(c(union_edges$from, union_edges$to)))

# Keep the comparison visually stable across reruns by anchoring the layout on
# the truth network rather than on the union of method-specific predictions.
layout_graph <- graph_from_data_frame(
  truth_edges[, .(from, to)],
  directed = FALSE,
  vertices = data.frame(name = vertices)
)
set.seed(42)
layout_xy <- layout_with_fr(layout_graph)
rownames(layout_xy) <- V(layout_graph)$name

plot_graph <- function(edge_dt, main_title, show_legend = FALSE) {
  g <- graph_from_data_frame(
    edge_dt,
    directed = FALSE,
    vertices = data.frame(name = vertices)
  )
  coords <- layout_xy[V(g)$name, , drop = FALSE]
  eweights <- E(g)$weight
  edge_widths <- if (length(eweights)) 1 + 5 * (eweights / max(eweights)) else numeric()
  edge_colors <- ifelse(
    E(g)$edge_sign == "positive", "#1f78b4",
    ifelse(E(g)$edge_sign == "negative", "#e31a1c", "grey55")
  )

  plot(
    g,
    layout = coords,
    vertex.label = V(g)$name,
    vertex.label.cex = 0.75,
    vertex.label.family = "sans",
    vertex.label.color = "black",
    vertex.size = 18,
    vertex.color = "grey96",
    vertex.frame.color = "grey35",
    vertex.label.dist = 0.2,
    edge.width = edge_widths,
    edge.color = edge_colors,
    edge.curved = 0.1,
    main = main_title,
    margin = 0.05
  )

  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], border = "grey45", lwd = 1.6, xpd = NA)

  if (show_legend) {
    legend(
      "bottomleft",
      legend = c("Positive", "Negative", "Mixed"),
      col = c("#1f78b4", "#e31a1c", "grey55"),
      lwd = c(2.5, 2.5, 2.5),
      bty = "n",
      cex = 0.78,
      inset = c(0.01, 0.02)
    )
  }
}

out_png <- file.path(fig_dir, sprintf("interaction_network_example_%s_%s.png", dataset, date_tag))
png(
  filename = out_png,
  width = 9.2,
  height = 3.8,
  units = "in",
  res = 300
)
op <- par(mfrow = c(1, 3), mar = c(0.3, 0.3, 2.1, 0.3), oma = c(0, 0, 0, 0))
plot_graph(truth_edges, "Truth", show_legend = TRUE)
plot_graph(pln_edges, "PLNNetwork")
plot_graph(glm_edges, "GLMNet (Poisson)")
par(op)
dev.off()

cat(out_png, "\n")
