library(data.table)
library(glmnet)
library(ggplot2)
library(MASS)

# Load and prepare data
dataname <- "glne007"
task.dt <- data.table::fread( paste0("~/Projects/pln_eval/data/", dataname, "_update.csv") )
#task.dt[, Group_ID := NULL]
task.dt <- log1p(task.dt)

# Find most abundant column and set as target
most_abundant_col <- names(which.max(colSums(task.dt, na.rm = TRUE)))
target_y <- task.dt[[most_abundant_col]]
print(paste("Target variable:", most_abundant_col))

# Remove the most abundant column from the feature matrix
task.dt[, (most_abundant_col) := NULL]

# Calculate correlations between each feature and the target
correlations <- task.dt[, lapply(.SD, function(x) cor(x, target_y, use = "complete.obs"))]
corr_df <- data.table(
  feature = names(correlations),
  correlation = as.numeric(correlations)
)
# Filter for positive correlations only
corr_df <- corr_df[correlation > 0]


# Use quantiles to create 5 equal-sized groups
corr_df[, corr_group := cut(correlation, 
                            breaks = quantile(correlation, probs = seq(0, 1, 0.2)),
                            labels = c("Lowest 20%", "Low 20%", "Mid 20%", "High 20%", "Highest 20%"),
                            include.lowest = TRUE)]

# Get the top 5 correlations from each group (should give exactly 25 features)
top5_corr <- corr_df[, head(.SD[order(-correlation)], 5), by = corr_group]

# Create correlation heatmap with better visibility for low correlations
p <- ggplot(top5_corr, aes(x = corr_group, y = reorder(feature, correlation), fill = correlation)) +
  geom_tile() +
  scale_fill_gradient(low = "lightgray", high = "red", name = "Correlation") +
  labs(title = paste("Top 5 Positive Correlations per Group with Target:", most_abundant_col),
       x = "Correlation Group", y = "Features") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

ggsave(paste0("~/Projects/pln_eval/out/", dataname, "_corr_heatmap.jpg"), p, width = 5.5, height = 4, dpi = 300)


# Create separate tasks for each correlation group
task.list <- list()
for (group in unique(top5_corr$corr_group)) {
  # Get features for this correlation group
  group_features <- top5_corr[corr_group == group]$feature
  # Create dataset with target and this group's features
  group_data <- data.table(
    target = target_y,
    task.dt[, ..group_features]
  )
  # Create mlr3 task for this group
  task_id <- paste0("corr_group_", gsub("[^A-Za-z0-9]", "_", group))
  reg_task <- mlr3::TaskRegr$new(
    id = task_id,
    backend = group_data,
    target = "target"
  )
  task.list[[task_id]] <- reg_task
  cat("Task:", task_id, "- Features:", reg_task$n_features, "\n")
}

