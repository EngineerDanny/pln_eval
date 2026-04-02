library(data.table)
library(mlr3)
library(atime)
library(PLNmodels)
library(glmnet)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
data_dir <- file.path(base_dir, "data")
out_dir <- file.path(base_dir, "out", "atime")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(base_dir, "load_source.R"))

if (length(args) < 1) {
  stop(
    "Usage: Rscript atime_pln_glmnet_p_scaling.R <dataname_or_path> ",
    "[n_fixed] [p_values_csv] [times] [seed] [pln_numepoch]"
  )
}

dataname_or_path <- args[1]
n_fixed <- if (length(args) >= 2) as.integer(args[2]) else 200L
p_values <- if (length(args) >= 3) {
  as.integer(strsplit(args[3], ",", fixed = TRUE)[[1]])
} else {
  c(10L, 25L, 50L, 100L, 150L, 200L)
}
times <- if (length(args) >= 4) as.integer(args[4]) else 3L
seed <- if (length(args) >= 5) as.integer(args[5]) else 1L
pln_numepoch <- if (length(args) >= 6) as.integer(args[6]) else 1000L

load_counts <- function(dataname_or_path, data_dir) {
  if (grepl("\\.tsv\\.gz$", dataname_or_path)) {
    tsv_path <- if (file.exists(dataname_or_path)) {
      dataname_or_path
    } else {
      file.path(data_dir, dataname_or_path)
    }
    if (!file.exists(tsv_path)) stop("File not found: ", dataname_or_path)
    raw <- read.table(
      gzfile(tsv_path),
      sep = "\t",
      header = TRUE,
      row.names = 1,
      check.names = FALSE
    )
    counts <- t(as.matrix(raw))
    counts <- counts[rowSums(counts) > 0, colSums(counts) > 0, drop = FALSE]
    return(as.data.table(counts))
  }

  tsv_path <- list.files(
    data_dir,
    pattern = paste0(dataname_or_path, ".*\\.tsv\\.gz$"),
    recursive = TRUE,
    full.names = TRUE
  )
  csv_path <- file.path(data_dir, paste0(dataname_or_path, "_update.csv"))

  if (length(tsv_path) > 0) {
    raw <- read.table(
      gzfile(tsv_path[1]),
      sep = "\t",
      header = TRUE,
      row.names = 1,
      check.names = FALSE
    )
    counts <- t(as.matrix(raw))
    counts <- counts[rowSums(counts) > 0, colSums(counts) > 0, drop = FALSE]
    as.data.table(counts)
  } else if (file.exists(csv_path)) {
    dt <- fread(csv_path)
    if ("Group_ID" %in% names(dt)) dt[, Group_ID := NULL]
    dt
  } else {
    stop("No data file found for: ", dataname_or_path)
  }
}

dataname <- if (grepl("\\.tsv\\.gz$", dataname_or_path)) {
  sub("\\.tsv\\.gz$", "", basename(dataname_or_path))
} else {
  dataname_or_path
}

counts_dt <- load_counts(dataname_or_path, data_dir)
if (nrow(counts_dt) < 2L || ncol(counts_dt) < 2L) {
  stop("Need at least 2 samples and 2 taxa after filtering.")
}

set.seed(seed)
n_use <- min(n_fixed, nrow(counts_dt))
row_idx <- sort(sample.int(nrow(counts_dt), n_use))
counts_dt <- counts_dt[row_idx]

taxa_rank <- names(sort(colSums(counts_dt), decreasing = TRUE))
max_p <- min(max(p_values), length(taxa_rank))
p_values <- sort(unique(p_values[p_values >= 2L & p_values <= max_p]))
if (length(p_values) < 2L) {
  stop("Need at least two valid p values after filtering.")
}

selected_taxa <- taxa_rank[seq_len(max_p)]
target_taxon <- selected_taxa[1]
target_label <- paste0("Taxa", target_taxon)

cat(
  sprintf(
    "Dataset=%s n=%d max_p=%d target=%s p_grid=%s times=%d\n",
    dataname,
    nrow(counts_dt),
    max_p,
    target_taxon,
    paste(p_values, collapse = ","),
    times
  )
)

time_result <- atime::atime(
  N = p_values,
  setup = {
    current_taxa <- selected_taxa[seq_len(N)]
    current_counts <- copy(counts_dt[, ..current_taxa])
    current_counts[, names(current_counts) := lapply(.SD, log1p)]
    setnames(current_counts, names(current_counts), paste0("Taxa", names(current_counts)))
    current_task <- mlr3::TaskRegr$new(target_label, current_counts, target = target_label)
  },
  expr.list = list(
    "GLMNet (Poisson)" = quote({
      x <- data.matrix(current_task$data(cols = current_task$feature_names))
      y <- current_task$data(cols = current_task$target_names)[[1]]
      glmnet::cv.glmnet(
        x,
        y,
        alpha = 1,
        family = "poisson",
        type.measure = "deviance"
      )
    }),
    "PLN" = quote({
      counts <- data.matrix(
        current_task$data(
          cols = c(current_task$target_names, current_task$feature_names)
        ),
        rownames.force = TRUE
      )
      rownames(counts) <- paste0("R", seq_len(nrow(counts)))
      torch::torch_set_num_threads(1L)
      PLN(
        Abundance ~ 1,
        data = pln_prepare_counts(counts, "none"),
        control = PLN_param(
          covariance = "full",
          trace = 0L,
          backend = "torch",
          config_optim = list(
            algorithm = "ADAM",
            lr = 0.01,
            numepoch = pln_numepoch
          )
        )
      )
    })
  ),
  times = times,
  seconds.limit = Inf,
  result = FALSE
)

meas_dt <- as.data.table(time_result$measurements)
meas_csv_dt <- meas_dt[, names(meas_dt)[!vapply(meas_dt, is.list, logical(1L))], with = FALSE]
fwrite(
  meas_csv_dt,
  file.path(out_dir, paste0(dataname, "_pln_glmnet_p_scaling_measurements.csv"))
)
saveRDS(
  time_result,
  file.path(out_dir, paste0(dataname, "_pln_glmnet_p_scaling_atime.rds"))
)

best_list <- atime::references_best(time_result)
best_dt <- as.data.table(best_list$meas)
best_csv_dt <- best_dt[, names(best_dt)[!vapply(best_dt, is.list, logical(1L))], with = FALSE]
fwrite(
  best_csv_dt,
  file.path(out_dir, paste0(dataname, "_pln_glmnet_p_scaling_best.csv"))
)

meas_dt <- as.data.table(time_result$measurements)
if (nrow(meas_dt) > 0) {
  time_plot_dt <- meas_dt[
    expr.name %in% c("GLMNet (Poisson)", "PLN"),
    .(
      unit = "seconds",
      expr.name,
      expr.label = fifelse(expr.name == "GLMNet (Poisson)", "GLMNet", expr.name),
      N,
      empirical = median,
      ymin = q25,
      ymax = q75
    )
  ]
  mem_plot_dt <- meas_dt[
    expr.name %in% c("GLMNet (Poisson)", "PLN"),
    .(
      unit = "megabytes",
      expr.name,
      expr.label = fifelse(expr.name == "GLMNet (Poisson)", "GLMNet", expr.name),
      N,
      empirical = kilobytes / 1024,
      ymin = NA_real_,
      ymax = NA_real_
    )
  ]
  plot_dt <- rbindlist(list(mem_plot_dt, time_plot_dt), use.names = TRUE)
  plot_dt[, unit := factor(unit, levels = c("seconds", "megabytes"))]
  unit_label_map <- c(
    seconds = "Wall time (s)",
    megabytes = "Peak memory (MB)"
  )
  scaling_plot <- ggplot(plot_dt, aes(N, empirical, color = expr.name, fill = expr.name)) +
    geom_ribbon(
      data = plot_dt[!is.na(ymin)],
      aes(ymin = ymin, ymax = ymax),
      alpha = 0.18,
      linewidth = 0
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    facet_wrap(
      ~ unit,
      scales = "free_y",
      ncol = 2,
      labeller = as_labeller(unit_label_map)
    ) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = c("GLMNet (Poisson)" = "#D55E00", "PLN" = "#1f78b4")) +
    scale_fill_manual(values = c("GLMNet (Poisson)" = "#D55E00", "PLN" = "#1f78b4")) +
    theme_bw() +
    labs(
      x = "p",
      y = "Median fit cost\nacross 3 runs\n(log scale)",
      color = "Method",
      fill = "Method"
    ) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(linewidth = 0.4),
      strip.background = element_rect(linewidth = 0.4)
    )
  if (requireNamespace("directlabels", quietly = TRUE)) {
    x_limits <- c(min(plot_dt$N), max(plot_dt$N) * 2.2)
    scaling_plot <- scaling_plot +
      directlabels::geom_dl(
        aes(label = expr.label),
        method = "right.polygons",
        data = plot_dt[expr.name == "PLN"],
        cex = 0.9
      ) +
      directlabels::geom_dl(
        aes(label = expr.label),
        method = "right.polygons",
        data = plot_dt[expr.name == "GLMNet (Poisson)"],
        cex = 1.2
      ) +
      coord_cartesian(xlim = x_limits, clip = "off") +
      theme(
        legend.position = "none",
        plot.margin = margin(5.5, 65, 5.5, 5.5)
      )
  }
  ggsave(
    filename = file.path(out_dir, paste0(dataname, "_pln_glmnet_p_scaling.png")),
    plot = scaling_plot,
    width = 8.5,
    height = 2.3,
    dpi = 300
  )
}

cat("Wrote:\n")
cat(file.path(out_dir, paste0(dataname, "_pln_glmnet_p_scaling_measurements.csv")), "\n")
cat(file.path(out_dir, paste0(dataname, "_pln_glmnet_p_scaling_best.csv")), "\n")
cat(file.path(out_dir, paste0(dataname, "_pln_glmnet_p_scaling_atime.rds")), "\n")
