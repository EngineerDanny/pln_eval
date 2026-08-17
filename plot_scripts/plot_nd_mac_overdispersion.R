suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
benchmark <- fread(file.path(
  base_dir, "out", "count_prediction_results_retained20_with_stats.csv"
))
plot_dt <- copy(benchmark)
plot_dt[, outcome := fifelse(winner == "PLN", "PLN wins", "GLMNet wins")]
plot_dt[, outcome := factor(
  outcome,
  levels = c("PLN wins", "GLMNet wins")
)]
p_scatter <- ggplot(
  plot_dt,
  aes(
    x = N_over_D,
    y = mac,
    color = outcome
  )
) +
  geom_point(shape = 16, size = 2.4, alpha = 0.9) +
  scale_x_log10(
    limits = c(0.05, 100),
    breaks = c(0.1, 0.3, 1, 3, 10, 30, 100)
  ) +
  scale_y_continuous(
    limits = c(0, 0.25),
    breaks = seq(0, 0.25, 0.05),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_color_manual(
    values = c("#0000BF", "#BF0000"),
    name = NULL
  ) +
  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(shape = 16, size = 3.2)
    )
  ) +
  labs(
    x = expression(N/D~"ratio (log scale)"),
    y = "Mean absolute correlation (MAC)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linetype = "dotted", color = "grey78"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.margin = margin(2, 0, 0, 0),
    legend.key.height = unit(0.9, "lines"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9, color = "black"),
    panel.border = element_rect(linewidth = 0.45)
  )

plot_dt[, relative_pln_advantage := (glmnet - pln) / glmnet]
plot_dt[, log_nd := log(N_over_D)]
plot_dt[, c("z_log_nd", "z_mac") := lapply(
  .SD, function(x) as.numeric(scale(x))
), .SDcols = c("log_nd", "mac")]

fit <- lm(
  relative_pln_advantage ~ z_log_nd + z_mac,
  data = plot_dt
)
X <- model.matrix(fit)
hc3_weights <- (residuals(fit) / (1 - hatvalues(fit)))^2
hc3_vcov <- solve(crossprod(X)) %*%
  crossprod(X, X * as.numeric(hc3_weights)) %*%
  solve(crossprod(X))
hc3_se <- sqrt(diag(hc3_vcov))
critical_value <- qt(0.975, df.residual(fit))
coefficient_dt <- data.table(
  term = names(coef(fit)),
  beta = coef(fit),
  lower = coef(fit) - critical_value * hc3_se,
  upper = coef(fit) + critical_value * hc3_se,
  p_value = 2 * pt(abs(coef(fit) / hc3_se), df.residual(fit), lower.tail = FALSE)
)[term != "(Intercept)"]
coefficient_dt[, predictor := c("log(N/D)", "MAC")]
fwrite(
  coefficient_dt,
  file.path(base_dir, "out", "dataset_characteristic_regression_hc3.csv")
)

table_dt <- coefficient_dt[, .(
  Predictor = predictor,
  `Beta` = sprintf("%.3f", beta),
  `95% CI` = sprintf("[%.3f, %.3f]", lower, upper),
  `p` = ifelse(
    p_value < 0.001,
    "<.001",
    sub("^0", "", sprintf("%.3f", p_value))
  )
)]
column_bounds <- c(0.02, 0.34, 0.52, 0.84, 0.98)
column_centers <- (head(column_bounds, -1) + tail(column_bounds, -1)) / 2
row_bounds <- c(0.50, 0.60, 0.70, 0.80)
row_centers <- (head(row_bounds, -1) + tail(row_bounds, -1)) / 2

p_table <- ggplot() +
  annotate(
    "rect", xmin = column_bounds[1], xmax = tail(column_bounds, 1),
    ymin = row_bounds[3], ymax = row_bounds[4],
    fill = "grey88", color = NA
  ) +
  annotate(
    "rect", xmin = column_bounds[1], xmax = tail(column_bounds, 1),
    ymin = row_bounds[1], ymax = row_bounds[2],
    fill = "grey96", color = NA
  ) +
  annotate(
    "segment",
    x = column_bounds, xend = column_bounds,
    y = min(row_bounds), yend = max(row_bounds),
    color = "grey45", linewidth = 0.35
  ) +
  annotate(
    "segment",
    x = min(column_bounds), xend = max(column_bounds),
    y = row_bounds, yend = row_bounds,
    color = "grey45", linewidth = 0.35
  ) +
  annotate(
    "text", x = column_centers, y = row_centers[3],
    label = names(table_dt), fontface = "bold", size = 2.55
  ) +
  annotate(
    "text", x = rep(column_centers, each = 2),
    y = rep(rev(row_centers[1:2]), times = 4),
    label = unlist(table_dt, use.names = FALSE), size = 2.2
  ) +
  annotate(
    "text", x = 0.03, y = 0.90, label = "B) Adjusted regression",
    hjust = 0, fontface = "bold", size = 3.2
  ) +
  annotate(
    "text", x = 0.03, y = 0.27,
    label = sprintf(
      "Outcome: relative PLN advantage\nRobust 95%% CI (n = 20)\nR\u00b2 = %.2f; adjusted R\u00b2 = %.2f",
      summary(fit)$r.squared, summary(fit)$adj.r.squared
    ),
    hjust = 0, vjust = 1, size = 2.55, lineheight = 1.1
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void()

p <- (p_scatter + labs(tag = "A)")) +
  p_table +
  plot_layout(widths = c(1.45, 1.05)) &
  theme(plot.tag = element_text(face = "bold", size = 11))

output_path <- file.path(
  base_dir, "paper", "sagmb", "figures", "figure2_nd_mac.png"
)
ggsave(output_path, p, width = 7.2, height = 3.6, dpi = 400, bg = "white")
cat("Saved:", output_path, "\n")
