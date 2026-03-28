args <- commandArgs(trailingOnly = TRUE)
dataset_spec <- if (length(args) >= 1) args[1] else "amgut2_update.csv"

read_dataset_spec <- function(dataset_spec) {
  if (length(dataset_spec) != 1L || !nzchar(dataset_spec)) {
    stop("dataset_spec must be a non-empty string")
  }

  if (file.exists(dataset_spec) && grepl("\\.(txt|list|lst)$", dataset_spec)) {
    lines <- trimws(readLines(dataset_spec, warn = FALSE))
    lines <- lines[nzchar(lines)]
    lines <- lines[!grepl("^#", lines)]
    if (length(lines) == 0L) stop("No datasets found in: ", dataset_spec)
    return(lines)
  }

  dataset_spec
}

dataset_values <- read_dataset_spec(dataset_spec)
n_folds <- 3L
n_methods <- 4L
cat(length(dataset_values) * n_folds * n_methods, "\n")
