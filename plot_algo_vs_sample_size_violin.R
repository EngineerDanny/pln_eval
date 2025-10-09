library(data.table)
library(ggplot2)

sample_sizes <- c(30, 50, 100, 200, 289)
plot_data_list <- list()

for (i in 1:length(sample_sizes)) {
  temp_data <- fread(paste0("/projects/genomic-ml/da2343/PLN/pln_eval/bmr/amgut1_", sample_sizes[i], "_09_11.csv"))
  temp_data[, n_samples := sample_sizes[i]]
  plot_data_list[[i]] <- temp_data
}

plot_data <- rbindlist(plot_data_list)

gg <- ggplot(plot_data, aes(x = factor(n_samples), y = deviance, fill = algorithm)) +
  geom_violin(position = position_dodge(width = 0.75), alpha = 0.4, trim = FALSE, width = 0.8, color = NA) +
  geom_boxplot(width = 0.25, position = position_dodge(width = 0.75),
               outlier.shape = NA, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point",
               shape = 21, size = 2.4, fill = "white", color = "black",
               position = position_dodge(width = 0.75), stroke = 0.9) +
  scale_fill_manual(values = c("GLMNet" = "#E69F00", "PLN" = "#009E73")) +
  labs(
    title = "GLMNet dominates with small n, PLN overtakes with large n (Amgut1)",
    x = "Sample Size",
    y = "Poisson deviance (↓ better)",
    fill = "Algorithm"
  ) +
  coord_cartesian(ylim = c(0, 4.5)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave("/projects/genomic-ml/da2343/PLN/pln_eval/out/amgut1_violin_boxplot_09_11.png",
       plot = gg,
       width = 8.5,
       height = 4.6,
       dpi = 300)
