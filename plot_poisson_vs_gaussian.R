data_name <- "mixmpln"
data_path <- paste0("/projects/genomic-ml/da2343/PLN/pln_eval/data/poisson_vs_gaussian/", data_name, ".csv")
output_path <- paste0("/projects/genomic-ml/da2343/PLN/pln_eval/out/poisson_vs_gaussian_", data_name, ".png")

plot_data <- data.table::fread(data_path)
plot_data_filtered <- plot_data[learner_id %in% c("regr.cv_glmnet", "regr.pln")]

plot_data_filtered[, algorithm := fcase(
  learner_id == "regr.cv_glmnet", "Poisson GLMNET",
  learner_id == "regr.pln", "Poisson Log Normal"
)]

score_summary <- plot_data_filtered[, .(
  mean_deviance = mean(regr.poisson_deviance),
  sd_deviance = 1.96 * sd(regr.poisson_deviance)/sqrt(.N),
  mean_rmse = mean(regr.rmse),
  sd_rmse = 1.96 * sd(regr.rmse)/sqrt(.N),
  n = .N
), by = .(algorithm)]

score_summary[, `:=`(
  lower_ci_deviance = mean_deviance - sd_deviance,
  upper_ci_deviance = mean_deviance + sd_deviance,
  lower_ci_rmse = mean_rmse - sd_rmse,
  upper_ci_rmse = mean_rmse + sd_rmse,
  deviance_label = sprintf("%.3f", mean_deviance),
  rmse_label = sprintf("%.3f", mean_rmse)
)]

library(ggplot2)

gg <- ggplot(score_summary, aes(x = mean_rmse, y = mean_deviance, color = algorithm)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower_ci_deviance, ymax = upper_ci_deviance), 
                width = 0.01, size = 1) +
  geom_errorbarh(aes(xmin = lower_ci_rmse, xmax = upper_ci_rmse), 
                 height = 0.01, size = 1) +
  geom_text(aes(label = algorithm), 
            vjust = -0.8, size = 3.5, fontface = "bold") +
  scale_color_manual(values = c("Poisson GLMNET" = "#e74c3c", "Poisson Log Normal" = "#3498db")) +
  labs(
    x = "RMSE (lower better)",
    y = "Mean Poisson Deviance (lower better)",
    title = "Trade-off: RMSE vs. Poisson Deviance",
    subtitle = "Poisson GLMNET vs. Poisson Log Normal with 95% Confidence Intervals",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(output_path,
       plot = gg,
       width = 5, 
       height = 4,
       dpi = 300)