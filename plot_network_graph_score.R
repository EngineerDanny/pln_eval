library(data.table)
library(ggplot2)

base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
in_dir <- file.path(base_dir, "out", "network_graph")
fig_dir <- file.path(base_dir, "figures", "march26")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

rep_path <- file.path(in_dir, "network_graph_pseudologlik_replicates.csv")
sum_path <- file.path(in_dir, "network_graph_pseudologlik_summary.csv")

rep_dt <- fread(rep_path)
sum_dt <- fread(sum_path)

label_map <- c(
  LOTO_PLN_glasso = "LOTO PLN+glasso",
  PLN_glasso = "PLN+glasso",
  LOTO_glmnet_CV1se = "LOTO glmnet",
  Baseline_diagonal = "Diagonal baseline"
)

order_levels <- sum_dt[order(-mean_loglik), method]
rep_dt[, method_label := factor(label_map[method], levels = label_map[order_levels])]

p <- ggplot(rep_dt, aes(x = loglik, y = method_label, color = method_label, fill = method_label)) +
  geom_boxplot(width = 0.55, alpha = 0.18, outlier.shape = NA, linewidth = 0.4) +
  geom_point(
    position = position_jitter(height = 0.08, width = 0),
    size = 1.9,
    alpha = 0.9
  ) +
  scale_color_manual(values = c(
    "LOTO PLN+glasso" = "#00BFC4",
    "PLN+glasso" = "#56B4E9",
    "LOTO glmnet" = "#D55E00",
    "Diagonal baseline" = "grey40"
  )) +
  scale_fill_manual(values = c(
    "LOTO PLN+glasso" = "#00BFC4",
    "PLN+glasso" = "#56B4E9",
    "LOTO glmnet" = "#D55E00",
    "Diagonal baseline" = "grey40"
  )) +
  theme_bw() +
  labs(
    x = "Held-out Gaussian pseudo-loglik",
    y = NULL
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.4),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9)
  )

out_path <- file.path(fig_dir, "network_graph_score_boxplot.png")
ggsave(out_path, p, width = 5.2, height = 2.8, dpi = 400)
cat(out_path, "\n")
