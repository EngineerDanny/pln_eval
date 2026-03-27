library(data.table)
library(ggplot2)

base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
figures_dir <- file.path(base_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(base_dir, "load_source.R"))

dataset_info <- data.table(
  dataset = c(
    "diabimmune_karelia_16s_family",
    "amgut1",
    "20190424.CosteaPI_2017.metaphlan_bugs_list.stool_genus",
    "mbqc_integrated_otus_family"
  ),
  dataset_name = c(
    "Diabimmune Karelia family",
    "amgut1",
    "CosteaPI 2017 stool",
    "MBQC integrated family"
  ),
  N = c(1584L, 289L, 279L, 18270L),
  D = c(55L, 127L, 179L, 235L)
)

read_case_study_data <- function(dataset) {
  load(file.path(base_dir, "bmr", paste0(dataset, ".RData")))
  score_dt <- as.data.table(bmr$score(poisson_measure))
  score_dt[, dataset := dataset]
  score_dt[, algorithm := fcase(
    learner_id == "regr.featureless", "Featureless",
    learner_id == "regr.lasso", "GLMNet (Poisson)",
    learner_id == "regr.pln", "PLN"
  )]
  score_dt <- score_dt[!is.na(algorithm)]

  score_summary <- score_dt[, .(
    mean_deviance = mean(regr.poisson_deviance),
    se_deviance = sd(regr.poisson_deviance) / sqrt(.N)
  ), by = .(dataset, algorithm)]

  pairwise_results <- score_dt[, .(
    featureless = regr.poisson_deviance[algorithm == "Featureless"],
    glmnet = regr.poisson_deviance[algorithm == "GLMNet (Poisson)"],
    pln = regr.poisson_deviance[algorithm == "PLN"]
  ), by = .(dataset, task_id, iteration)]

  glmnet_test <- pairwise_results[, {
    test <- t.test(pln, glmnet, paired = TRUE)
    .(
      mean_diff = mean(glmnet - pln),
      p_value = test$p.value
    )
  }, by = dataset]

  rm(bmr, score_dt, pairwise_results)
  gc()

  list(score_summary = score_summary, glmnet_test = glmnet_test)
}

dataset_results <- lapply(dataset_info$dataset, read_case_study_data)
score_summary <- rbindlist(lapply(dataset_results, `[[`, "score_summary"))

score_summary[, `:=`(
  lower = mean_deviance - se_deviance,
  upper = mean_deviance + se_deviance,
  label = sprintf("%.2f±%.02f", mean_deviance, se_deviance)
)]

score_summary <- merge(score_summary, dataset_info, by = "dataset", all.x = TRUE)
score_summary[, dataset_label := paste0(dataset_name, "\nn=", N, ", p=", D)]
score_summary[, dataset_label := factor(
  dataset_label,
  levels = paste0(dataset_info$dataset_name, "\nn=", dataset_info$N, ", p=", dataset_info$D)
)]

score_summary[, algorithm := factor(algorithm, levels = c("Featureless", "GLMNet (Poisson)", "PLN"))]

glmnet_tests <- rbindlist(lapply(dataset_results, `[[`, "glmnet_test"))

glmnet_tests[, p_label := ifelse(
  p_value < 1e-4,
  "P<1e-4",
  sprintf("P=%.4f", p_value)
)]
glmnet_tests[, diff_label := sprintf("Diff=%.3f, %s", mean_diff, p_label)]

glmnet_ref <- score_summary[algorithm == "GLMNet (Poisson)", .(
  dataset,
  glmnet_mean = mean_deviance
)]

pln_annotations <- merge(
  score_summary[algorithm == "PLN", .(dataset, algorithm, mean_deviance, lower, dataset_label)],
  glmnet_ref,
  by = "dataset",
  all.x = TRUE
)
pln_annotations <- merge(pln_annotations, glmnet_tests, by = "dataset", all.x = TRUE)

segment_data <- pln_annotations[, .(
  dataset,
  dataset_label,
  algorithm,
  mean_deviance,
  lower,
  glmnet_mean,
  diff_label
)]

gg <- ggplot(score_summary, aes(x = mean_deviance, y = algorithm)) +
  geom_segment(
    data = segment_data,
    aes(x = glmnet_mean, xend = mean_deviance, y = algorithm, yend = algorithm),
    inherit.aes = FALSE,
    color = "#D55E00",
    alpha = 0.6,
    linewidth = 1.2
  ) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    orientation = "y",
    width = 0.12,
    linewidth = 0.3
  ) +
  geom_point(shape = 1, size = 0.7) +
  geom_text(
    aes(label = label),
    color = "black",
    vjust = 1.5,
    size = 1.5
  ) +
  geom_text(
    data = segment_data,
    aes(x = lower, y = algorithm, label = diff_label),
    inherit.aes = FALSE,
    color = "#D55E00",
    hjust = 0,
    vjust = -0.65,
    size = 1.5
  ) +
  facet_wrap(~ dataset_label, nrow = 1, scales = "free_x") +
  scale_y_discrete(limits = c("Featureless", "GLMNet (Poisson)", "PLN")) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  labs(
    x = "Mean Poisson Deviance (Mean ± SE)",
    y = "Algorithm"
  ) +
  coord_cartesian(clip = "off") +
  theme_grey() +
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

out_file <- file.path(figures_dir, "glmnet_case_study_benchmark_4datasets.png")
ggsave(out_file, gg, width = 6.7, height = 1.7, dpi = 700)
cat("Saved:", out_file, "\n")
