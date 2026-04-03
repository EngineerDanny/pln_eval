## Plot F1 bar chart with bootstrap SE error bars for Figure 4
## Requires: out/bootstrap_f1_results.csv (from bootstrap_f1.R)

library(data.table)
library(ggplot2)

base_dir  <- "/projects/genomic-ml/da2343/PLN/pln_eval"
fig_dir   <- file.path(base_dir, "paper", "archive")

bootstrap_path <- file.path(base_dir, "out", "bootstrap_f1_results.csv")
if (!file.exists(bootstrap_path)) stop("Run bootstrap_f1.R first.")

results <- fread(bootstrap_path)

dataset_info <- data.table(
  dataset = c(
    "omm12",
    "omm12_keystone_2023",
    "pairinterax",
    "butyrate_assembly_2021",
    "host_fitness_2018"
  ),
  dataset_label = c(
    "OMM12",
    "OMM12\nkeystone 2023",
    "PairInteraX",
    "Butyrate\nassembly 2021",
    "Host fitness\n2018"
  ),
  benchmark_group = c(
    "Broad edge recovery",
    "Broad edge recovery",
    "Broad edge recovery",
    "Local / directional",
    "Local / directional"
  )
)

# Compute mean and SE per dataset × method from bootstrap
summary_dt <- results[, .(
  mean_f1 = mean(f1, na.rm = TRUE),
  se_f1   = sd(f1, na.rm = TRUE) / sqrt(sum(!is.na(f1)))
), keyby = .(dataset, method)]

# Merge in labels
plot_dt <- merge(summary_dt, dataset_info, by = "dataset", all.x = TRUE)
plot_dt <- plot_dt[method %in% c("PLNnetwork", "LOTO_glmnet_CV1se")]

method_map <- c(PLNnetwork = "PLNNetwork", LOTO_glmnet_CV1se = "GLMNet (Poisson)")
plot_dt[, method_label := factor(method_map[method], levels = c("PLNNetwork", "GLMNet (Poisson)"))]
plot_dt[, benchmark_group := factor(benchmark_group, levels = c("Broad edge recovery", "Local / directional"))]
plot_dt[, dataset_label := factor(dataset_label, levels = dataset_info$dataset_label)]

fill_map <- c("PLNNetwork" = "grey30", "GLMNet (Poisson)" = "grey80")

p <- ggplot(plot_dt, aes(x = dataset_label, y = mean_f1, fill = method_label,
                         ymin = pmax(mean_f1 - se_f1, 0), ymax = pmin(mean_f1 + se_f1, 1))) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.66,
    color = "black",
    linewidth = 0.3
  ) +
  geom_errorbar(
    position = position_dodge(width = 0.75),
    width = 0.25,
    linewidth = 0.5
  ) +
  facet_grid(. ~ benchmark_group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = fill_map, name = NULL) +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(x = NULL, y = "Edge-recovery F1") +
  theme_classic(base_size = 10) +
  theme(
    legend.position    = "top",
    legend.text        = element_text(size = 10),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 10, face = "plain"),
    axis.text.x        = element_text(size = 9, color = "black"),
    axis.text.y        = element_text(size = 9, color = "black"),
    axis.title.y       = element_text(size = 11, margin = margin(r = 6)),
    axis.line          = element_line(linewidth = 0.4),
    axis.ticks         = element_line(linewidth = 0.4),
    panel.spacing.x    = unit(1.1, "lines"),
    plot.margin        = margin(6, 8, 6, 6)
  )

out_png <- file.path(fig_dir, "figure4_network_f1_benchmarks.png")
ggsave(out_png, p, width = 7.1, height = 3.2, dpi = 400)
cat("Saved:", out_png, "\n")

# Print summary for reference
cat("\n=== Bootstrap F1 mean ± SE ===\n")
print(plot_dt[, .(dataset, method_label, mean_f1 = round(mean_f1, 3), se_f1 = round(se_f1, 3))][order(dataset, method_label)])
