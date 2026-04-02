library(data.table)
library(igraph)

args <- commandArgs(trailingOnly = TRUE)
export_tag <- if (length(args) >= 1) args[1] else "top_20"
top_edges_per_method <- if (length(args) >= 2) as.integer(args[2]) else 25L

base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
export_dir <- file.path(base_dir, "out", "network_graph", export_tag)
fig_dir <- file.path(base_dir, "figures", "march26")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

edge_path <- file.path(export_dir, "network_full_fit_edges.csv")
node_path <- file.path(export_dir, "network_full_fit_nodes.csv")

if (!file.exists(edge_path)) stop("Missing edge file: ", edge_path)
if (!file.exists(node_path)) stop("Missing node file: ", node_path)

edge_dt <- fread(edge_path)
node_dt <- fread(node_path)

method_levels <- c("LOTO_PLN_glasso", "PLN_glasso", "LOTO_glmnet_CV1se")
method_labels <- c(
  LOTO_PLN_glasso = "LOTO PLN+glasso",
  PLN_glasso = "PLN+glasso",
  LOTO_glmnet_CV1se = "LOTO glmnet"
)

edge_dt <- edge_dt[method %in% method_levels]
edge_top <- edge_dt[, head(.SD[order(-abs_weight)], top_edges_per_method), by = method]

union_pairs <- unique(edge_top[, .(taxon_i, taxon_j)])
vertices <- unique(c(union_pairs$taxon_i, union_pairs$taxon_j))
node_plot <- node_dt[taxon %in% vertices]
node_plot <- node_plot[match(vertices, taxon)]

union_graph <- graph_from_data_frame(union_pairs, directed = FALSE, vertices = node_plot[, .(name = taxon)])
set.seed(42)
layout_xy <- layout_with_fr(union_graph)
rownames(layout_xy) <- V(union_graph)$name

plot_method_graph <- function(method_name) {
  sub_edges <- edge_top[method == method_name]
  g <- graph_from_data_frame(
    sub_edges[, .(from = taxon_i, to = taxon_j, weight, abs_weight, sign)],
    directed = FALSE,
    vertices = node_plot[, .(
      name = taxon,
      total_abundance,
      prevalence,
      mean_log1p,
      abundance_rank
    )]
  )

  coords <- layout_xy[V(g)$name, , drop = FALSE]
  vertex_sizes <- 5 + 10 * (V(g)$total_abundance / max(V(g)$total_abundance))
  edge_widths <- 1 + 5 * (E(g)$abs_weight / max(E(g)$abs_weight))
  edge_cols <- ifelse(E(g)$sign == "positive", "#00BFC4", "#D55E00")

  plot(
    g,
    layout = coords,
    vertex.label = V(g)$name,
    vertex.label.cex = 0.7,
    vertex.label.family = "sans",
    vertex.label.color = "black",
    vertex.size = vertex_sizes,
    vertex.color = "grey95",
    vertex.frame.color = "grey30",
    vertex.label.dist = 0.35,
    edge.width = edge_widths,
    edge.color = edge_cols,
    edge.curved = 0.08,
    main = method_labels[[method_name]],
    margin = 0.05
  )
}

png(
  filename = file.path(fig_dir, sprintf("representative_network_%s.png", export_tag)),
  width = 10,
  height = 3.3,
  units = "in",
  res = 300
)
op <- par(mfrow = c(1, 3), mar = c(0.4, 0.4, 1.8, 0.4), oma = c(0, 0, 0, 0))
for (method_name in method_levels) {
  plot_method_graph(method_name)
}
par(op)
dev.off()

cat(file.path(fig_dir, sprintf("representative_network_%s.png", export_tag)), "\n")
