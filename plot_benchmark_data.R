
#datanames <- c("amgut2_297_samples", "ioral_86_samples", "baxter_crc_490_samples") # PLN better
#datanames <- c("crohns_100_samples", "amgut1_289_samples", "mixmpln_195_samples") # Glmnet better
#datanames <- c("hmp2prot_47_samples", "hmp216S_47_samples")

library(data.table)
library(ggplot2)
library(ggrepel)

dataset_info <- data.table(
  dataset = c("hmp2prot_47_samples", "hmp216S_47_samples"),
  dataset_name = c("hmp2prot", "hmp216S"),
  N = c(47, 47),
  D = c(43, 45)
)

dataset_info <- data.table(
  dataset = c( "hmp216S_47_samples"),
  dataset_name = c("hmp216S"),
  N = c(47),
  D = c( 45)
)

#dataset_info <- data.table(
#  dataset = c("crohns_100_samples", "amgut1_289_samples", "mixmpln_195_samples"),
#  dataset_name = c("crohns", "amgut1", "mixmpln"),
#  N = c(100, 289, 195),
#  D = c(5, 127, 129)
#)

#dataset_info <- data.table(
#  dataset = c("amgut2_297_samples", "ioral_86_samples", "baxter_crc_490_samples"),
#  dataset_name = c("amgut2", "ioral", "baxter_crc"),
#  N = c(296, 86, 490),
#  D = c(138, 63, 117)
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
  sd_deviance = 1.96 * sd(regr.poisson_deviance)/sqrt(.N),
  n = .N
), by = .(algorithm, dataset)]
score_summary[, `:=`(
  lower_std = mean_deviance - sd_deviance,
  upper_std = mean_deviance + sd_deviance,
  label = sprintf("%.2f", mean_deviance)
)]
score_summary <- merge(score_summary, dataset_info, by = "dataset")
score_summary$dataset <- factor(score_summary$dataset, levels = dataset_info$dataset)
score_summary[, dataset_clean := gsub("_06_20|_06_22", "", dataset)]
score_summary[, dataset_label := paste0("data: ", dataset_name, "\nN=", N, ", D=", D)]
score_summary$dataset_clean <- factor(score_summary$dataset_clean, levels = gsub("_06_20|_06_22", "", dataset_info$dataset))
algorithm_order <- c("Featureless", "GLMNet (Poisson)", "PLN", "PLN Network", "PLN PCA")
score_summary$algorithm <- factor(score_summary$algorithm, levels = algorithm_order)
score_summary[, is_pln := algorithm %in% c("PLN", "PLN Network", "PLN PCA")]
gg <- ggplot(score_summary, aes(x = mean_deviance, y = algorithm)) +
  geom_errorbarh(aes(xmin = lower_std, xmax = upper_std), 
                 height = 0.1, linewidth = 0.3) +
  geom_point( shape = 1, size = 0.7) +
  geom_text(aes(label = label), 
            nudge_x = 0.2, 
            hjust = 0,
            size = 1.5) +
  facet_wrap(~ dataset_label, nrow = 2, ncol = 3, scales = "free_x") +
  labs(
    title = "PLN PCA outperforms others",
    #title = "GLMNet (Poisson) marginally outperforms PLN",
    #title = "PLN models outperform GLMNet (Poisson)",
    x = "Mean Poisson Deviance (95% CI)",
    y = "Algorithm"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  theme(
    axis.text.y = element_text(color = "black", size = 5),
    axis.text.x = element_text(size = 5),
    strip.text = element_text(size = 6),
    axis.title = element_text(size = 6),
    plot.title = element_text(size = 7),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
ggsave(
       "/projects/genomic-ml/da2343/PLN/pln_eval/data/poisson_vs_gaussian/PLN_special_bmr.png",
       #"/projects/genomic-ml/da2343/PLN/pln_eval/data/poisson_vs_gaussian/glmnet_outperforms_bmr.png",
       plot = gg,
       width = 3, #4, 
       height = 1.4, #1.4,
       dpi = 300)
