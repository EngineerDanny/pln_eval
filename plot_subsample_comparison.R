library(data.table)
library(ggplot2)

base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
source(file.path(base_dir, "load_source.R"))

make_subsample_plot <- function(dataset_prefix, dataset_label, sample_sizes, bmr_dir, out_dir, log_scale = FALSE) {
  plot_data_list <- list()
  for (n in sample_sizes) {
    rdata_file <- file.path(bmr_dir, paste0(dataset_prefix, "_", n, "_samples.RData"))
    if (!file.exists(rdata_file)) {
      cat("Skipping (not found):", rdata_file, "\n")
      next
    }
    cat("Loading:", rdata_file, "\n")
    load(rdata_file)  # loads 'bmr'
    score.dt <- bmr$score(poisson_measure)
    # Keep only lasso and PLN (drop featureless)
    score.dt <- score.dt[learner_id %in% c("regr.lasso", "regr.pln")]
    score.dt[, algorithm := ifelse(learner_id == "regr.lasso", "GLMNet", "PLN")]
    score.dt[, n_samples := n]
    plot_data_list[[length(plot_data_list) + 1]] <- score.dt[, .(n_samples, algorithm, deviance = regr.poisson_deviance)]
  }

  if (length(plot_data_list) == 0) {
    cat("No data found for", dataset_label, "\n")
    return(NULL)
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
    labs(
      title = paste0(dataset_label, ": GLMNet vs PLN by Sample Size"),
      x = "Sample Size",
      y = if (log_scale) "log10(Poisson deviance)" else "Poisson deviance (Lower is better)",
      fill = "Algorithm"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )

  if (log_scale) {
    gg <- gg + scale_y_log10()
  }

  suffix <- if (log_scale) "_log" else ""
  out_file <- file.path(out_dir, paste0(dataset_label, "_algo_vs_sample_size", suffix, ".png"))
  ggsave(out_file, plot = gg, width = 8.5, height = 4.6, dpi = 300)
  cat("Saved:", out_file, "\n")
  gg
}

bmr_dir <- file.path(base_dir, "bmr")
out_dir <- file.path(base_dir, "out")

make_subsample_plot("deblur_125nt_no_blooms_family", "Deblur_family",
                    c(100, 500, 1000, 2000), bmr_dir, out_dir)

make_subsample_plot("mbqc_integrated_otus_family", "MBQC_family",
                    c(100, 500, 1000, 2000), bmr_dir, out_dir)

make_subsample_plot("emp_cr_gg_13_8.release1_family", "EMP_family",
                    c(100, 500, 1000, 2000), bmr_dir, out_dir)

make_subsample_plot("amgut2", "amgut2",
                    c(30, 50, 100, 200, 297), bmr_dir, out_dir)

## Log-scale versions
make_subsample_plot("deblur_125nt_no_blooms_family", "Deblur_family",
                    c(100, 500, 1000, 2000), bmr_dir, out_dir, log_scale = TRUE)

make_subsample_plot("mbqc_integrated_otus_family", "MBQC_family",
                    c(100, 500, 1000, 2000), bmr_dir, out_dir, log_scale = TRUE)

make_subsample_plot("emp_cr_gg_13_8.release1_family", "EMP_family",
                    c(100, 500, 1000, 2000), bmr_dir, out_dir, log_scale = TRUE)

make_subsample_plot("amgut2", "amgut2",
                    c(30, 50, 100, 200, 297), bmr_dir, out_dir, log_scale = TRUE)
