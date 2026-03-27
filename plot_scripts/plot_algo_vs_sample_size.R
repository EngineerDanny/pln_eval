library(data.table)
library(ggplot2)

figures_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval/figures"
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

sample_sizes <- c(30, 50, 100, 200, 297)
plot_data_list <- list()

for (i in 1:length(sample_sizes)) {
  temp_data <- fread(paste0("/projects/genomic-ml/da2343/PLN/pln_eval/bmr/amgut2_", sample_sizes[i], "_09_11.csv"))
  temp_data[, n_samples := sample_sizes[i]]
  plot_data_list[[i]] <- temp_data
}

plot_data <- rbindlist(plot_data_list)

gg <- ggplot(plot_data, aes(x = factor(n_samples), y = deviance, fill = algorithm)) +
  geom_boxplot(width = 0.55, alpha = 0.55, outlier.shape = NA,
               position = position_dodge2(width = 0.6)) +
  geom_point(aes(color = algorithm),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
             alpha = 0.25, size = 0.9, shape = 16) +
  scale_fill_manual(values = c("GLMNet" = "#E69F00", "PLN" = "#009E73")) +
  scale_color_manual(values = c("GLMNet" = "#E69F00", "PLN" = "#009E73"), guide = "none") +
  coord_cartesian(ylim = c(0, 4.5)) +
  labs(
    title = "GLMNet dominates with small N, PLN overtakes with large N (Amgut2)",
    x = "Sample Size",
    y = "Poisson deviance (Lower is better)",
    fill = "Algorithm"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(figures_dir, "amgut2_boxplot_jittered_refined_09_11.png"),
       plot = gg,
       width = 8.5,
       height = 4.6,
       dpi = 300)
