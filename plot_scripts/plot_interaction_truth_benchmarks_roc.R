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

read_tsv_gz <- function(path) {
  read.delim(gzfile(path), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

sort_pair <- function(a, b) {
  left <- ifelse(a <= b, a, b)
  right <- ifelse(a <= b, b, a)
  data.table(taxon_1 = left, taxon_2 = right)
}

roc_curve_dt <- function(score, truth) {
  pos_total <- sum(truth == 1L)
  neg_total <- sum(truth == 0L)

  if (pos_total == 0L || neg_total == 0L) {
    return(data.table(fpr = c(0, 1), tpr = c(0, 1)))
  }

  thresholds <- sort(unique(score), decreasing = TRUE)
  thresholds <- c(Inf, thresholds, -Inf)

  rows <- lapply(thresholds, function(thr) {
    pred_pos <- score >= thr
    tp <- sum(pred_pos & truth == 1L)
    fp <- sum(pred_pos & truth == 0L)
    data.table(
      fpr = fp / neg_total,
      tpr = tp / pos_total
    )
  })

  unique(rbindlist(rows))
}

auc_trapz <- function(fpr, tpr) {
  ord <- order(fpr, tpr)
  x <- fpr[ord]
  y <- tpr[ord]
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
  LOTO_glmnet_CV1se = "LOTO glmnet"
)

roc_rows <- list()
auc_rows <- list()
k <- 1L
for (i in seq_len(nrow(dataset_info))) {
  dataset <- dataset_info$dataset[i]
  for (method in method_levels) {
    scored <- build_scored_pairs(dataset, method)
    roc_dt <- roc_curve_dt(scored$score, scored$truth)
    roc_dt[, `:=`(
      dataset = dataset,
      method = method,
      auc = auc_trapz(fpr, tpr)
    )]
    roc_rows[[k]] <- roc_dt
    auc_rows[[k]] <- unique(roc_dt[, .(dataset, method, auc)])
    k <- k + 1L
  }
}

plot_dt <- rbindlist(roc_rows)
auc_dt <- rbindlist(auc_rows)
plot_dt <- merge(plot_dt, dataset_info, by = "dataset", all.x = TRUE)
plot_dt[, method_label := factor(
  method_map[method],
  levels = c("PLNNetwork", "LOTO glmnet")
)]
plot_dt[, dataset_label := factor(dataset_label, levels = dataset_info$dataset_label)]

fill_map <- c(
  "PLNNetwork" = "grey20",
  "LOTO glmnet" = "grey75"
)
line_map <- c(
  "PLNNetwork" = "solid",
  "LOTO glmnet" = "dotdash"
)

p <- ggplot(plot_dt, aes(x = fpr, y = tpr, color = method_label, linetype = method_label)) +
  geom_abline(intercept = 0, slope = 1, color = "grey85", linewidth = 0.35) +
  geom_line(linewidth = 0.55) +
  facet_wrap(~ dataset_label, nrow = 2) +
  scale_color_manual(values = fill_map, name = NULL) +
  scale_linetype_manual(values = line_map, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.5), expand = c(0, 0)) +
  labs(
    x = "False positive rate",
    y = "True positive rate"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 10),
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 11),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    plot.margin = margin(6, 8, 6, 6)
  )

out_png <- file.path(fig_dir, "interaction_truth_benchmarks_roc.png")
ggsave(out_png, p, width = 7.1, height = 4.2, dpi = 400)

fwrite(
  merge(auc_dt, dataset_info[, .(dataset, dataset_label, benchmark_group)], by = "dataset", all.x = TRUE)[
    , .(dataset, dataset_label, benchmark_group, method, auc)
  ],
  file.path(fig_dir, "interaction_truth_benchmarks_roc_auc.tsv"),
  sep = "\t"
)

cat(out_png, "\n")
