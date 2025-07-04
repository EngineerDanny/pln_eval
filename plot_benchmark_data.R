library(data.table)
library(ggplot2)


datanames <- c( "amgut2_06_20", "ioral_06_22","baxter_crc_06_22", 
                "amgut1_06_20", "hmpv35_06_20", "mixmpln_06_22")

datanames <- c( "crohns_06_22", "amgut1_06_20", "hmpv35_06_20")
plot_data_list <- list()
for(i in 1:length(datanames)) {
  temp_data <- fread(paste0("/projects/genomic-ml/da2343/PLN/pln_eval/data/plot_data/", datanames[i], ".csv"))
  temp_data[, dataset := datanames[i]]
  plot_data_list[[i]] <- temp_data
}
plot_data <- rbindlist(plot_data_list)
plot_data[, algorithm := fcase(
  learner_id == "regr.featureless", "Featureless",
  learner_id == "regr.cv_glmnet", "GLMNet (Poisson)",
  learner_id == "regr.pln", "PLN"
)]
score_summary <- plot_data[, .(
  mean_deviance = mean(regr.poisson_deviance),
  sd_deviance = 1.96 * sd(regr.poisson_deviance)/sqrt(.N),
  n = .N
), by = .(algorithm, dataset)]
score_summary[, `:=`(
  lower_ci = mean_deviance - sd_deviance,
  upper_ci = mean_deviance + sd_deviance,
  label = sprintf("%.3f", mean_deviance)
)]
p_values <- data.table()
for(ds in datanames) {
  baseline_values <- plot_data[learner_id == "regr.featureless" & dataset == ds, regr.poisson_deviance]
  for (alg in c("regr.cv_glmnet", "regr.pln")) {
    alg_values <- plot_data[learner_id == alg & dataset == ds, regr.poisson_deviance]
    
    if (length(alg_values) == length(baseline_values)) {
      t_result <- t.test(alg_values, baseline_values, paired = TRUE)
    } else {
      t_result <- t.test(alg_values, baseline_values, paired = FALSE)
    }
    
    clean_alg <- ifelse(alg == "regr.cv_glmnet", "GLMNet (Poisson)", "PLN")
    
    p_values <- rbind(p_values, data.table(
      algorithm = clean_alg,
      dataset = ds,
      p_value = t_result$p.value
    ))
  }
  p_values <- rbind(p_values, data.table(
    algorithm = "Featureless",
    dataset = ds,
    p_value = 1.0
  ))
}
score_summary <- merge(score_summary, p_values, by = c("algorithm", "dataset"))
score_summary$dataset <- factor(score_summary$dataset, levels = datanames)
score_summary[, dataset_clean := gsub("_06_20|_06_22", "", dataset)]
score_summary$dataset_clean <- factor(score_summary$dataset_clean, levels = gsub("_06_20|_06_22", "", datanames))
score_summary[, p_value_label := fcase(
  p_value >= 1.0, "",
  p_value < 0.001, "p < 0.001",
  p_value < 0.01, "p < 0.01", 
  p_value < 0.05, "p < 0.05",
  default = paste0("p = ", round(p_value, 3))
)]
algorithm_order <- c("Featureless", "GLMNet (Poisson)", "PLN")
score_summary$algorithm <- factor(score_summary$algorithm, levels = algorithm_order)
score_summary[, is_pln := algorithm == "PLN"]
gg <- ggplot(score_summary, aes(x = mean_deviance, y = algorithm)) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, color = is_pln), 
                 height = 0.1, linewidth = 0.5) +
  geom_point(aes(color = is_pln), size = 1) +
  geom_text(aes(label = label, color = is_pln), 
            vjust = 1.8, size = 3.5) +
  geom_text(data = score_summary[p_value_label != ""], 
            aes(label = p_value_label, color = is_pln), 
            vjust = -1.1, size = 3, fontface = "italic") +
  scale_color_manual(values = c("black", "red"), 
                     labels = c("Other Methods", "PLN")) +
  facet_wrap(~ dataset_clean, nrow = 2, ncol = 3) +
  labs(
    #title = "PLN Shows Consistent Superior Performance Across Datasets",
    title = "GLMNet marginally outperforms PLN",
    x = "Mean Poisson Deviance (95% CI)",
    y = "Algorithm"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(
      color = ifelse(levels(score_summary$algorithm) == "PLN", "red", "black")
    ),
    axis.text.x = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
ggsave("/projects/genomic-ml/da2343/PLN/pln_eval/out/2_bmr.png",
       plot = gg,
       width = 8, 
       #height = 4,
       height = 3,
       dpi = 300)