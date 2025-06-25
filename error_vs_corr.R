library(ggplot2)
library(data.table)

dataname <- "mixmpln"
# Read the results
results_data <- fread( paste0( "~/Projects/pln_eval/data/", dataname, "_corr_benchmark.csv" ) )

# Filter for PLN and cv_glmnet only (exclude featureless baseline)
plot_data <- results_data[learner_id %in% c("regr.pln", "regr.cv_glmnet")]

# Fix the correlation group extraction
plot_data[, corr_group := gsub("corr_group_", "", task_id)]
plot_data[, corr_group := gsub("_", " ", corr_group)]
plot_data[, corr_group := gsub(" $", "%", corr_group)]  # Replace trailing space with %

# Create ordered factor for correlation groups
plot_data[, corr_group := factor(corr_group,
                                 levels = c("Lowest 20%", "Low 20%", "Mid 20%", "High 20%", "Highest 20%"))]

# Map learner names using data.table syntax
plot_data[learner_id == "regr.pln", method := "Poisson Log Normal (Covariance-Aware)"]
plot_data[learner_id == "regr.cv_glmnet", method := "Poisson GLMNET (Independent)"]


# Create the plot with all individual points
library(directlabels)


p <- ggplot(plot_data, aes(x = corr_group, 
                           y = regr.poisson_deviance, color = method)) +
  geom_point(size = 2) +
  geom_line(aes(group = method), stat = "summary", fun = mean, size = 1.2) +
  geom_dl(aes(label = method), method = list("smart.grid", cex = 0.6)) +
  scale_color_manual(values = c("Poisson Log Normal (Covariance-Aware)" = "#1f77b4", 
                                "Poisson GLMNET (Independent)" = "#ff7f0e")) +
  labs(x = "Correlation Groups",
       y = "Poisson deviance error",
       title = dataname,
       color = "") +
  theme(legend.position = "none")


ggsave(paste0("~/Projects/pln_eval/out/", dataname, "_error_vs_corr.jpg"), p, width = 4.5, height = 3.5, dpi = 300)

