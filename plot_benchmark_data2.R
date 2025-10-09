library(data.table)
library(ggplot2)
library(ggrepel)

dataset_info <- data.table(
  dataset = c("amgut2_297_samples", "ioral_86_samples", "baxter_crc_490_samples"),
  dataset_name = c("amgut2", "ioral", "baxter_crc"),
  N = c(296, 86, 490),
  D = c(138, 63, 117)
)

dataset_info <- data.table(
  dataset = c( "hmp216S_47_samples"),
  dataset_name = c("hmp216S"),
  N = c(47),
  D = c(45)
)

#dataset_info <- data.table(
#  dataset = c("crohns_100_samples", "mixmpln_195_samples"),
#  dataset_name = c("crohns", "mixmpln"),
#  N = c(100, 195),
#  D = c(5, 129)
#)


plot_data_list <- list()
for(i in 1:nrow(dataset_info)) {
  temp_data <- fread(paste0("/projects/genomic-ml/da2343/PLN/pln_eval/data/poisson_vs_gaussian/", dataset_info$dataset[i], ".csv"))
  temp_data[, dataset := dataset_info$dataset[i]]
  plot_data_list[[i]] <- temp_data
}

plot_data <- rbindlist(plot_data_list)

plot_data[, algorithm := fcase(
  learner_id == "regr.featureless", "Featureless",
  learner_id == "regr.cv_glmnet", "GLMNet (Poisson)",
  learner_id == "regr.pln", "PLN",
  learner_id == "regr.pln_network.tuned", "PLN Network",
  learner_id == "regr.plnpca.tuned", "PLN PCA"
)]

score_summary <- plot_data[, .(
  mean_deviance = mean(regr.poisson_deviance),
  se_deviance = sd(regr.poisson_deviance) / sqrt(.N),
  n = .N
), by = .(algorithm, dataset)]

score_summary[, `:=`(
  lower_std = mean_deviance - se_deviance,
  upper_std = mean_deviance + se_deviance,
  label = sprintf("%.2f±%.2f", mean_deviance, se_deviance)
)]

score_summary <- merge(score_summary, dataset_info, by = "dataset")
score_summary$dataset <- factor(score_summary$dataset, levels = dataset_info$dataset)
score_summary[, dataset_clean := gsub("_06_20|_06_22", "", dataset)]
score_summary[, dataset_label := paste0("data: ", dataset_name, "\nN=", N, ", D=", D)]

dataset_label_levels <- score_summary[order(match(dataset, dataset_info$dataset)), unique(dataset_label)]
score_summary$dataset_label <- factor(score_summary$dataset_label, levels = dataset_label_levels)

score_summary$dataset_clean <- factor(score_summary$dataset_clean, levels = gsub("_06_20|_06_22", "", dataset_info$dataset))

algorithm_order <- c("Featureless", "GLMNet (Poisson)", "PLN", "PLN Network", "PLN PCA")
score_summary$algorithm <- factor(score_summary$algorithm, levels = algorithm_order)
score_summary[, is_pln := algorithm %in% c("PLN", "PLN Network", "PLN PCA")]

baseline <- "GLMNet (Poisson)"

baseline_data <- plot_data[algorithm == baseline]
setnames(baseline_data, "regr.poisson_deviance", "baseline_deviance")
baseline_data <- baseline_data[, .(dataset, baseline_deviance)]

comparison_data <- merge(plot_data, baseline_data, by = "dataset", allow.cartesian = TRUE)

test_results <- comparison_data[algorithm != baseline, {
  test <- t.test(regr.poisson_deviance, baseline_deviance, paired = TRUE)
  browser()
  .(mean_diff = mean(regr.poisson_deviance - baseline_deviance),
    p_value = test$p.value)
}, by = .(algorithm, dataset)]

baseline_means <- score_summary[algorithm == baseline, .(dataset, baseline_mean = mean_deviance)]

score_summary_with_tests <- merge(score_summary, test_results, by = c("algorithm", "dataset"), all.x = TRUE)
score_summary_with_tests <- merge(score_summary_with_tests, baseline_means, by = "dataset", all.x = TRUE)

segment_data <- score_summary_with_tests[!algorithm %in% c(baseline, "Featureless")]
segment_data[, performance := ifelse(mean_diff < 0, "Better", "Worse")]
segment_data[, p_label := ifelse(p_value < 0.0001, "P<0.0001", sprintf("P=%.4f", p_value))]
segment_data[, label_hjust := ifelse(mean_diff < 0, 0, 1)]

label_position <- segment_data[, .(algorithm, dataset, label_hjust)]
score_summary_with_tests <- merge(score_summary_with_tests, label_position, by = c("algorithm", "dataset"), all.x = TRUE)
score_summary_with_tests[is.na(label_hjust), label_hjust := 0.5]

gg <- ggplot(score_summary_with_tests, aes(x = mean_deviance, y = algorithm)) +
  geom_segment(aes(x = baseline_mean, xend = mean_deviance, 
                   y = algorithm, yend = algorithm, color = performance),
               alpha = 0.5, linewidth = 1.5,
               data = segment_data) +
  geom_text(aes(x = mean_deviance, y = algorithm,
                label = sprintf("Diff=%.4f, %s", mean_diff, p_label),
                hjust = label_hjust,
                color = performance),
            size = 1.5, vjust = -0.5,
            data = segment_data) +
  scale_color_manual(values = c("Better" = "#0072B2", "Worse" = "#D55E00")) +
  geom_errorbarh(aes(xmin = lower_std, xmax = upper_std), 
                 height = 0.1, linewidth = 0.3) +
  geom_point(shape = 1, size = 0.7) +
  geom_text(aes(label = label, hjust = label_hjust), 
            vjust = 1.5,
            size = 1.5) +
  facet_wrap(~ dataset_label, nrow = 2, ncol = 3, scales = "free_x") +
  coord_cartesian(clip = "off") +
  labs(
    x = "Mean Poisson Deviance (Mean ± SE)",
    y = "Algorithm"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  theme(
    axis.text.y = element_text(color = "black", size = 5),
    axis.text.x = element_text(size = 5),
    strip.text = element_text(size = 6),
    axis.title = element_text(size = 6),
    plot.title = element_text(size = 7),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(1, 1, 2, 1), "lines")
  )

ggsave(
  "/projects/genomic-ml/da2343/PLN/pln_eval/data/poisson_vs_gaussian/plnpca_wins.png",
  plot = gg,
  #width = 5,
  #height = 2,
  width = 2.5,
  height = 1.8,
  dpi = 700)
