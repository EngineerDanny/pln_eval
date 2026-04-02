library(data.table)
library(ggplot2)

base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
truth_dir <- file.path(base_dir, "data", "interaction_ground_truth")
fig_dir <- file.path(base_dir, "figures", "march26")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
date_tag <- format(Sys.Date(), "%Y-%m-%d")

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
selected_dataset <- "omm12_keystone_2023"

read_tsv_gz <- function(path) {
  read.delim(gzfile(path), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

sort_pair <- function(a, b) {
  left <- ifelse(a <= b, a, b)
  right <- ifelse(a <= b, b, a)
  data.table(taxon_1 = left, taxon_2 = right)
}

pr_curve_dt <- function(score, truth) {
  pos_total <- sum(truth == 1L)
  if (pos_total == 0L) {
    return(data.table(recall = c(0, 1), precision = c(1, 1)))
  }

  ord <- order(score, decreasing = TRUE)
  score <- score[ord]
  truth <- truth[ord]

  tp <- cumsum(truth == 1L)
  fp <- cumsum(truth == 0L)
  keep <- c(diff(score) != 0, TRUE)

  curve <- data.table(
    recall = tp[keep] / pos_total,
    precision = tp[keep] / (tp[keep] + fp[keep])
  )

  curve <- rbind(
    data.table(recall = 0, precision = 1),
    curve,
    data.table(
      recall = 1,
      precision = mean(truth)
    )
  )
  unique(curve[order(recall, -precision)])
}

auprc_trapz <- function(recall, precision) {
  ord <- order(recall, precision)
  x <- recall[ord]
  y <- precision[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

build_scored_pairs <- function(dataset, method) {
  proc_dir <- file.path(truth_dir, dataset, "processed")
  truth <- as.data.table(read_tsv_gz(file.path(proc_dir, "truth_undirected.tsv.gz")))
  if (dataset %in% c("pairinterax", "butyrate_assembly_2021", "host_fitness_2018")) {
    abundance <- as.data.table(read_tsv_gz(file.path(proc_dir, "abundance_matrix.tsv.gz")))
  } else if (dataset %in% c("omm12", "omm12_keystone_2023")) {
    abundance <- as.data.table(read_tsv_gz(file.path(proc_dir, "community_absabundance_in_vitro.tsv.gz")))
  } else {
    stop("Unsupported dataset: ", dataset)
  }
  taxa <- setdiff(names(abundance), "sample_id")

  all_pairs_mat <- t(combn(sort(taxa), 2))
  all_pairs <- data.table(
    taxon_1 = all_pairs_mat[, 1],
    taxon_2 = all_pairs_mat[, 2]
  )
  all_pairs[, key := paste(taxon_1, taxon_2, sep = "||")]

  truth_pairs <- copy(truth[, .(taxon_1, taxon_2)])
  truth_pairs[, key := paste(taxon_1, taxon_2, sep = "||")]
  truth_keys <- unique(truth_pairs$key)

  pred_path <- file.path(proc_dir, "benchmark_outputs", paste0(method, "_predictions.tsv.gz"))
  pred <- as.data.table(read_tsv_gz(pred_path))
  if (nrow(pred)) {
    pred_pairs <- cbind(sort_pair(pred$taxon_1, pred$taxon_2), score = abs(pred$weight))
    pred_pairs[, key := paste(taxon_1, taxon_2, sep = "||")]
    pred_pairs <- pred_pairs[, .(score = max(score, na.rm = TRUE)), by = key]
  } else {
    pred_pairs <- data.table(key = character(), score = numeric())
  }

  scored <- merge(all_pairs[, .(key)], pred_pairs, by = "key", all.x = TRUE)
  scored[is.na(score), score := 0]
  scored[, truth := as.integer(key %in% truth_keys)]
  scored
}

method_levels <- c("PLNnetwork", "LOTO_glmnet_CV1se")
method_map <- c(
  PLNnetwork = "PLNNetwork",
  LOTO_glmnet_CV1se = "GLMNet (Poisson)"
)

curve_rows <- list()
auc_rows <- list()
k <- 1L
dataset <- selected_dataset
for (method in method_levels) {
  scored <- build_scored_pairs(dataset, method)
  pr_dt <- pr_curve_dt(scored$score, scored$truth)
  prevalence <- mean(scored$truth)
  pr_dt[, `:=`(
    dataset = dataset,
    method = method,
    auprc = auprc_trapz(recall, precision),
    prevalence = prevalence
  )]
  curve_rows[[k]] <- pr_dt
  auc_rows[[k]] <- unique(pr_dt[, .(dataset, method, auprc, prevalence)])
  k <- k + 1L
}

plot_dt <- rbindlist(curve_rows)
auc_dt <- rbindlist(auc_rows)
plot_dt <- merge(plot_dt, dataset_info, by = "dataset", all.x = TRUE)
plot_dt[, method_label := factor(
  method_map[method],
  levels = c("PLNNetwork", "GLMNet (Poisson)")
)]
plot_dt[, dataset_label := factor(dataset_label, levels = dataset_info[dataset == selected_dataset]$dataset_label)]

baseline_dt <- unique(plot_dt[, .(dataset_label, prevalence)])

line_map <- c(
  "PLNNetwork" = "solid",
  "GLMNet (Poisson)" = "dotdash"
)

p <- ggplot(plot_dt, aes(x = recall, y = precision, color = method_label, linetype = method_label)) +
  geom_hline(
    data = baseline_dt,
    aes(yintercept = prevalence),
    inherit.aes = FALSE,
    color = "grey85",
    linewidth = 0.35
  ) +
  geom_line(linewidth = 0.55) +
  scale_color_discrete(name = NULL) +
  scale_linetype_manual(values = line_map, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.5), expand = c(0, 0)) +
  labs(
    x = "Recall",
    y = "Precision",
    title = dataset_info[dataset == selected_dataset]$dataset_label
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 11),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    plot.margin = margin(6, 8, 6, 6)
  )

out_png <- file.path(fig_dir, sprintf("interaction_truth_benchmarks_auprc_%s_%s.png", selected_dataset, date_tag))
ggsave(out_png, p, width = 4.3, height = 3.5, dpi = 400)

fwrite(
  merge(auc_dt, dataset_info[, .(dataset, dataset_label, benchmark_group)], by = "dataset", all.x = TRUE)[
    , .(dataset, dataset_label, benchmark_group, method, auprc, prevalence)
  ],
  file.path(fig_dir, sprintf("interaction_truth_benchmarks_auprc_%s_%s.tsv", selected_dataset, date_tag)),
  sep = "\t"
)

cat(out_png, "\n")
