library(data.table)
library(ggplot2)

base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
truth_dir <- file.path(base_dir, "data", "interaction_ground_truth")
fig_dir <- file.path(base_dir, "figures", "march26")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

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

read_metrics <- function(dataset) {
  path <- file.path(truth_dir, dataset, "processed", "benchmark_outputs", "metrics_summary.tsv")
  dt <- fread(path)
  dt[, dataset := dataset]
  dt
}

plot_dt <- rbindlist(lapply(dataset_info$dataset, read_metrics), fill = TRUE)
plot_dt <- merge(plot_dt, dataset_info, by = "dataset", all.x = TRUE)
plot_dt <- plot_dt[method %in% c("PLNnetwork", "LOTO_glmnet_CV1se")]

method_map <- c(
  PLNnetwork = "PLNNetwork",
  LOTO_glmnet_CV1se = "LOTO glmnet"
)
plot_dt[, method_label := factor(
  method_map[method],
  levels = c("PLNNetwork", "LOTO glmnet")
)]
plot_dt[, benchmark_group := factor(
  benchmark_group,
  levels = c("Broad edge recovery", "Local / directional")
)]

dataset_levels <- dataset_info$dataset_label
plot_dt[, dataset_label := factor(dataset_label, levels = dataset_levels)]

fill_map <- c(
  "PLNNetwork" = "grey30",
  "LOTO glmnet" = "grey80"
)

p <- ggplot(plot_dt, aes(x = dataset_label, y = f1, fill = method_label)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.66,
    color = "black",
    linewidth = 0.3
  ) +
  facet_grid(. ~ benchmark_group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = fill_map, name = NULL) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = NULL,
    y = "Edge-recovery F1"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 10),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "plain"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.y = element_text(size = 11, margin = margin(r = 6)),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    panel.spacing.x = unit(1.1, "lines"),
    plot.margin = margin(6, 8, 6, 6)
  )

out_png <- file.path(fig_dir, "interaction_truth_benchmarks_f1.png")
ggsave(out_png, p, width = 7.1, height = 3.0, dpi = 400)

cat(out_png, "\n")
