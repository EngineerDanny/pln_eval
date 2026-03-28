library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
dataset_tag <- if (length(args) >= 1) args[1] else "amgut2_update"
in_dir <- file.path(base_dir, "out", "network_graph", dataset_tag)
fig_dir <- file.path(base_dir, "figures", "march26")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

rep_path <- file.path(in_dir, "network_graph_pseudologlik_replicates.csv")
rep_dt <- fread(rep_path)

label_map <- c(
  LOTO_PLN_glasso   = "LOTO PLN+glasso",
  PLN_glasso        = "PLN+glasso",
  LOTO_glmnet_CV1se = "LOTO glmnet",
  Baseline_diagonal = "Diagonal baseline"
)

baseline_dt <- rep_dt[method == "Baseline_diagonal", .(rep, baseline_loglik = loglik)]
plot_dt <- merge(rep_dt, baseline_dt, by = "rep", all.x = TRUE)
plot_dt[, improvement_pct := 100 * (loglik - baseline_loglik) / abs(baseline_loglik)]

sum_dt <- plot_dt[, .(
  mean_improvement_pct = mean(improvement_pct),
  sd_improvement_pct = sd(improvement_pct)
), by = method]
sum_dt <- sum_dt[method != "Baseline_diagonal"]

# Best model at top after coord_flip
order_levels <- sum_dt[order(mean_improvement_pct), method]
sum_dt[, method_label := factor(label_map[method], levels = label_map[order_levels])]

p <- ggplot(sum_dt, aes(x = method_label, y = mean_improvement_pct)) +
  geom_col(
    width = 0.62,
    fill = "grey70",
    color = "black",
    linewidth = 0.35
  ) +
  geom_errorbar(
    aes(
      ymin = mean_improvement_pct - sd_improvement_pct,
      ymax = mean_improvement_pct + sd_improvement_pct
    ),
    width = 0.2,
    linewidth = 0.5,
    color = "black"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed", color = "grey40") +
  coord_flip() +
  theme_classic(base_size = 10) +
  labs(
    x = NULL,
    y = "Improvement over diagonal baseline (%)"
  ) +
  theme(
    axis.text          = element_text(size = 10, color = "black"),
    axis.title.x       = element_text(size = 11, margin = margin(t = 6)),
    axis.line          = element_line(linewidth = 0.4),
    axis.ticks         = element_line(linewidth = 0.4),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    plot.margin        = margin(6, 10, 6, 6)
  )

# PNG for inspection
out_png <- file.path(fig_dir, sprintf("network_graph_score_improvement_pct_%s.png", dataset_tag))
ggsave(out_png, p, width = 5.2, height = 2.8, dpi = 400)

# PDF for submission
out_pdf <- file.path(fig_dir, sprintf("network_graph_score_improvement_pct_%s.pdf", dataset_tag))
ggsave(out_pdf, p, width = 5.2, height = 2.8)

cat(out_png, "\n")
