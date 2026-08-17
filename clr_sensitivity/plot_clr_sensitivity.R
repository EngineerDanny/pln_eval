library(data.table)
library(ggplot2)

base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
input_file <- file.path(base_dir, "clr_sensitivity", "case_study_summary.dat")
output_file <- file.path(
  base_dir, "paper", "sagmb", "figures", "figure_clr_sensitivity.png"
)

plot_dt <- fread(input_file)
test_dt <- fread(file.path(base_dir, "clr_sensitivity", "paired_tests.dat"))
plot_dt[, method := factor(
  method,
  levels = c("Featureless", "CLR-GLMNet", "GLMNet")
)]
plot_dt[, panel := dataset_label]
plot_dt[, panel := factor(panel, levels = unique(panel))]
plot_dt[, `:=`(
  lower = mean_deviance - se_deviance,
  upper = mean_deviance + se_deviance,
  label = sprintf("%.2f\u00b1%.02f", mean_deviance, se_deviance)
)]
panel_limits <- plot_dt[, .(
  panel_min = min(lower),
  panel_max = max(upper)
), by = panel]
plot_dt <- merge(plot_dt, panel_limits, by = "panel", sort = FALSE)
plot_dt[, `:=`(
  label_x = ifelse(mean_deviance > (panel_min + panel_max) / 2, upper, lower),
  label_hjust = ifelse(mean_deviance > (panel_min + panel_max) / 2, 1, 0)
)]

test_dt[, diff_label := ifelse(
  p_value < 1e-4,
  sprintf("Diff=%.3f, p<1e-4", diff_clr_minus_glmnet),
  sprintf("Diff=%.3f, p=%.4f", diff_clr_minus_glmnet, p_value)
)]
test_dt <- merge(
  test_dt,
  unique(plot_dt[, .(dataset, panel)]),
  by = "dataset"
)
method_means <- dcast(
  plot_dt[method %in% c("GLMNet", "CLR-GLMNet")],
  dataset + panel ~ method,
  value.var = "mean_deviance"
)
test_dt <- merge(test_dt, method_means, by = c("dataset", "panel"))
test_dt[, panel_color := "#0B559F"]

p <- ggplot(plot_dt, aes(x = mean_deviance, y = method)) +
  geom_segment(
    data = test_dt,
    aes(x = GLMNet, xend = `CLR-GLMNet`, y = "CLR-GLMNet", yend = "CLR-GLMNet"),
    inherit.aes = FALSE,
    color = test_dt$panel_color,
    alpha = 0.55,
    linewidth = 1.4
  ) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    orientation = "y",
    width = 0.12,
    linewidth = 0.35
  ) +
  geom_point(shape = 1, size = 2.2, stroke = 0.7) +
  geom_text(
    aes(x = label_x, label = label, hjust = label_hjust),
    vjust = 1.5,
    size = 3
  ) +
  geom_text(
    data = test_dt,
    aes(x = pmin(GLMNet, `CLR-GLMNet`), y = "CLR-GLMNet", label = diff_label),
    inherit.aes = FALSE,
    color = test_dt$panel_color,
    hjust = 0,
    vjust = -0.8,
    size = 3
  ) +
  facet_wrap(~panel, nrow = 1, scales = "free_x") +
  scale_y_discrete(
    limits = c("Featureless", "CLR-GLMNet", "GLMNet")
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.06, 0.10))
  ) +
  labs(
    x = "Held-out Poisson deviance (mean \u00b1 SE)",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", linewidth = 0.4),
    strip.text = element_text(face = "plain", size = 9),
    panel.border = element_rect(linewidth = 0.4)
  )

ggsave(output_file, p, width = 7.2, height = 2.2, dpi = 300, bg = "white")
cat("Saved:", output_file, "\n")
