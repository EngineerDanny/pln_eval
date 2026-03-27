library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
dataset_tag <- if (length(args) >= 1) args[1] else "amgut2_update"
in_dir <- file.path(base_dir, "out", "network_graph", dataset_tag)
fig_dir <- file.path(base_dir, "figures", "march26")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

rep_path <- file.path(in_dir, "network_graph_pseudologlik_replicates.csv")
sum_path <- file.path(in_dir, "network_graph_pseudologlik_summary.csv")

sum_dt <- fread(sum_path)

label_map <- c(
  LOTO_PLN_glasso = "LOTO PLN+glasso",
  PLN_glasso = "PLN+glasso",
  LOTO_glmnet_CV1se = "LOTO glmnet",
  Baseline_diagonal = "Diagonal baseline"
)

order_levels <- sum_dt[order(-mean_loglik), method]
sum_dt[, method_label := factor(label_map[method], levels = label_map[order_levels])]

p <- ggplot(sum_dt, aes(x = method_label, y = mean_loglik)) +
  geom_col(width = 0.62, fill = "grey75", color = "black", linewidth = 0.35) +
  geom_errorbar(
    aes(ymin = mean_loglik - sd_loglik, ymax = mean_loglik + sd_loglik),
    width = 0.16,
    linewidth = 0.35,
    color = "black"
  ) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Held-out Gaussian pseudo-loglik"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.4),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9)
  )

out_path <- file.path(fig_dir, sprintf("network_graph_score_boxplot_%s.png", dataset_tag))
ggsave(out_path, p, width = 5.2, height = 2.8, dpi = 400)
cat(out_path, "\n")
