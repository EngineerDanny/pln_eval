library(data.table)
library(ggplot2)
library(grid)

base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
figures_dir <- file.path(base_dir, "figures", "march26")
plot_data_file <- file.path(base_dir, "data", "processed", "count_prediction_case_study_benchmark_data.csv")
build_data_only <- "--build-data-only" %in% commandArgs(trailingOnly = TRUE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(plot_data_file), recursive = TRUE, showWarnings = FALSE)

source(file.path(base_dir, "load_source.R"))

dataset_info <- data.table(
  dataset = c(
    "amgut2",
    "crc_zeller",
    "hiv_lozupone_family",
    "cdi_schubert_family",
    "diabimmune_karelia_16s_family",
    "amgut1",
    "20190424.CosteaPI_2017.metaphlan_bugs_list.stool_genus",
    "mbqc_integrated_otus_family"
  ),
  dataset_name = c(
    "amgut2",
    "Zeller CRC genus",
    "Lozupone HIV family",
    "Schubert CDI family",
    "Diabimmune Karelia family",
    "amgut1",
    "CosteaPI 2017 stool",
    "MBQC integrated family"
  ),
  N = c(296L, 128L, 56L, 346L, 1584L, 289L, 279L, 18270L),
  D = c(138L, 257L, 62L, 79L, 55L, 127L, 179L, 235L),
  panel = c(
    "Representative datasets where PLN outperforms GLMNet (Poisson)",
    "Representative datasets where PLN outperforms GLMNet (Poisson)",
    "Representative datasets where PLN outperforms GLMNet (Poisson)",
    "Representative datasets where PLN outperforms GLMNet (Poisson)",
    "Representative datasets where GLMNet (Poisson) outperforms PLN",
    "Representative datasets where GLMNet (Poisson) outperforms PLN",
    "Representative datasets where GLMNet (Poisson) outperforms PLN",
    "Representative datasets where GLMNet (Poisson) outperforms PLN"
  ),
  panel_order = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
  panel_color = c(
    "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9",
    "#D55E00", "#D55E00", "#D55E00", "#D55E00"
  )
)

read_case_study_data <- function(dataset) {
  message("Loading ", dataset)
  data_env <- new.env(parent = emptyenv())
  load(file.path(base_dir, "bmr", paste0(dataset, ".RData")), envir = data_env)
  score_dt <- as.data.table(data_env$bmr$score(poisson_measure))
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
    glmnet = regr.poisson_deviance[algorithm == "GLMNet (Poisson)"],
    pln = regr.poisson_deviance[algorithm == "PLN"]
  ), by = .(dataset, task_id, iteration)]

  test_dt <- pairwise_results[, {
    test <- t.test(pln, glmnet, paired = TRUE)
    .(
      diff_pln_minus_glmnet = mean(pln - glmnet),
      diff_glmnet_minus_pln = mean(glmnet - pln),
      p_value = test$p.value
    )
  }, by = dataset]

  rm(data_env, score_dt, pairwise_results)
  gc()

  list(score_summary = score_summary, test_dt = test_dt)
}

build_plot_data <- function() {
  score_summary_list <- vector("list", length = nrow(dataset_info))
  test_dt_list <- vector("list", length = nrow(dataset_info))

  for (i in seq_len(nrow(dataset_info))) {
    dataset_result <- read_case_study_data(dataset_info$dataset[[i]])
    score_summary_list[[i]] <- dataset_result$score_summary
    test_dt_list[[i]] <- dataset_result$test_dt
    rm(dataset_result)
    gc()
  }

  score_summary <- rbindlist(score_summary_list)

  score_summary[, `:=`(
    lower = mean_deviance - se_deviance,
    upper = mean_deviance + se_deviance,
    label = sprintf("%.2f±%.02f", mean_deviance, se_deviance)
  )]

  score_summary <- merge(score_summary, dataset_info, by = "dataset", all.x = TRUE)
  score_summary[, dataset_label := paste0(dataset_name, "\nn=", N, ", p=", D)]

  panel_tests <- merge(
    unique(dataset_info[, .(dataset, panel, panel_color, panel_order)]),
    rbindlist(test_dt_list),
    by = "dataset",
    all.x = TRUE
  )
  rm(score_summary_list, test_dt_list)
  gc()

  panel_tests[, diff_value := ifelse(panel_order == 1L, diff_pln_minus_glmnet, diff_glmnet_minus_pln)]
  panel_tests[, p_label := ifelse(p_value < 1e-4, "P<1e-4", sprintf("P=%.4f", p_value))]
  panel_tests[, diff_label := sprintf("Diff=%.3f, %s", diff_value, p_label)]

  glmnet_ref <- score_summary[algorithm == "GLMNet (Poisson)", .(
    dataset,
    glmnet_mean = mean_deviance
  )]

  segment_data <- merge(
    score_summary[algorithm == "PLN", .(
      dataset, panel, panel_order, panel_color, mean_deviance, lower, dataset_label
    )],
    glmnet_ref,
    by = "dataset",
    all.x = TRUE
  )
  segment_data <- merge(
    segment_data,
    panel_tests[, .(dataset, diff_value, diff_label)],
    by = "dataset",
    all.x = TRUE
  )

  plot_dt <- merge(
    score_summary,
    segment_data[, .(
      dataset,
      glmnet_mean,
      diff_value,
      diff_label
    )],
    by = "dataset",
    all.x = TRUE
  )
  plot_dt[]
}

if (file.exists(plot_data_file)) {
  message("Using cached plot data from ", plot_data_file)
  plot_dt <- fread(plot_data_file)
} else {
  message("Building plot data and caching to ", plot_data_file)
  plot_dt <- build_plot_data()
  fwrite(plot_dt, plot_data_file)
}

if (build_data_only) {
  cat("Saved plot data:", plot_data_file, "\n")
  quit(save = "no", status = 0)
}

plot_dt[, algorithm := factor(algorithm, levels = c("Featureless", "GLMNet (Poisson)", "PLN"))]

build_panel_plot <- function(panel_name, show_y_title = TRUE) {
  panel_levels <- dataset_info[panel == panel_name][order(dataset_name), paste0(dataset_name, "\nn=", N, ", p=", D)]

  panel_scores <- copy(plot_dt[panel == panel_name])
  panel_scores[, dataset_label := factor(dataset_label, levels = panel_levels)]

  panel_segments <- unique(copy(plot_dt[panel == panel_name & algorithm == "PLN", .(
    dataset,
    panel,
    panel_order,
    panel_color,
    mean_deviance,
    lower,
    dataset_label,
    glmnet_mean,
    diff_label
  )]))
  panel_segments[, dataset_label := factor(dataset_label, levels = panel_levels)]

  panel_color <- unique(panel_scores$panel_color)

  ggplot(panel_scores, aes(x = mean_deviance, y = algorithm)) +
    geom_segment(
      data = panel_segments,
      aes(x = glmnet_mean, xend = mean_deviance, y = algorithm, yend = algorithm),
      inherit.aes = FALSE,
      color = panel_color,
      alpha = 0.55,
      linewidth = 1.4
    ) +
    geom_errorbar(
      aes(xmin = lower, xmax = upper),
      orientation = "y",
      width = 0.12,
      linewidth = 0.35
    ) +
    geom_point(shape = 1, size = 1.8, stroke = 0.6) +
    geom_text(
      aes(label = label),
      color = "black",
      vjust = 1.55,
      size = 2.3
    ) +
    geom_text(
      data = panel_segments,
      aes(x = lower, y = algorithm, label = diff_label),
      inherit.aes = FALSE,
      color = panel_color,
      size = 2.4,
      hjust = 0,
      vjust = -0.65
    ) +
    facet_wrap(~ dataset_label, nrow = 1, scales = "free_x") +
    scale_y_discrete(limits = c("Featureless", "GLMNet (Poisson)", "PLN")) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.08))) +
    labs(
      title = panel_name,
      x = "Mean Poisson Deviance (Mean ± SE)",
      y = if (show_y_title) "Algorithm" else NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_grey() +
    theme(
      axis.text.y = element_text(color = "black", size = 10),
      axis.text.x = element_text(size = 10),
      strip.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 12, hjust = 0, face = "plain"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(0.7, 0.8, 0.9, 0.8), "lines")
    )
}

top_plot <- build_panel_plot(
  "Representative datasets where PLN outperforms GLMNet (Poisson)",
  show_y_title = TRUE
)
bottom_plot <- build_panel_plot(
  "Representative datasets where GLMNet (Poisson) outperforms PLN",
  show_y_title = TRUE
)
out_file <- file.path(figures_dir, "count_prediction_case_study_benchmark.png")
message("Saving combined figure to ", out_file)
gc()
png(filename = out_file, width = 3600, height = 1700, res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(
  nrow = 3,
  ncol = 1,
  heights = unit(c(0.48, 0.48, 0.04), "npc")
)))
print(top_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
gc()
print(bottom_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.text(
  "Mean Poisson Deviance (Mean \u00b1 SE)",
  x = unit(0.5, "npc"),
  y = unit(0.5, "npc"),
  gp = gpar(fontsize = 14),
  vp = viewport(layout.pos.row = 3, layout.pos.col = 1)
)
dev.off()
cat("Saved:", out_file, "\n")
