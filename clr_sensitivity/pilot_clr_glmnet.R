library(glmnet)

base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"

load_dataset <- function(name) {
  if (name == "amgut1") {
    dat <- read.csv(file.path(base_dir, "data/amgut1_update.csv"), check.names = FALSE)
    if ("Group_ID" %in% names(dat)) dat$Group_ID <- NULL
    return(as.matrix(dat))
  }

  if (name == "hiv_lozupone_family") {
    path <- file.path(
      base_dir,
      "data/family/microbiomehd/hiv_lozupone_family.tsv.gz"
    )
    raw <- read.table(
      gzfile(path), sep = "\t", header = TRUE, row.names = 1,
      check.names = FALSE
    )
    counts <- t(as.matrix(raw))
    return(counts[rowSums(counts) > 0, colSums(counts) > 0, drop = FALSE])
  }

  stop("Unknown pilot dataset: ", name)
}

poisson_deviance <- function(truth, response) {
  response <- pmax(response, 1e-10)
  log_term <- ifelse(truth == 0, 0, truth * log(truth / response))
  mean(2 * (log_term - (truth - response)))
}

fit_dataset <- function(name, outer_folds = 3L, seed = 1L) {
  counts <- load_dataset(name)
  storage.mode(counts) <- "double"
  z <- log1p(counts)

  set.seed(seed)
  fold_id <- sample(rep(seq_len(outer_folds), length.out = nrow(z)))
  target_indices <- head(order(apply(z, 2, var), decreasing = TRUE), 5L)
  results <- vector("list", length(target_indices) * outer_folds * 2L)
  result_id <- 1L

  for (target in target_indices) {
    y <- z[, target]
    x_log1p <- z[, -target, drop = FALSE]
    x_clr <- x_log1p - rowMeans(x_log1p)

    for (fold in seq_len(outer_folds)) {
      train <- fold_id != fold
      test <- !train

      for (method in c("log1p", "target_excluded_clr")) {
        x <- if (method == "log1p") x_log1p else x_clr
        target_mean <- mean(y[train])
        fit <- tryCatch(
          suppressWarnings(cv.glmnet(
            x[train, , drop = FALSE], y[train], family = "poisson",
            alpha = 1, type.measure = "deviance", nfolds = 3
          )),
          error = function(e) NULL
        )
        pred <- if (is.null(fit)) {
          rep(target_mean, sum(test))
        } else {
          as.numeric(predict(
            fit, x[test, , drop = FALSE], s = "lambda.min", type = "response"
          ))
        }
        pred[!is.finite(pred) | pred < 0] <- target_mean
        pred <- pmin(pred, max(y[train]) * 2)

        results[[result_id]] <- data.frame(
          dataset = name,
          target = target,
          target_name = colnames(z)[target],
          fold = fold,
          method = method,
          deviance = poisson_deviance(y[test], pred)
        )
        result_id <- result_id + 1L
      }
    }
  }

  do.call(rbind, results)
}

datasets <- "hiv_lozupone_family"
detail <- do.call(rbind, lapply(datasets, fit_dataset))
summary <- aggregate(deviance ~ dataset + target + target_name + method, detail, mean)
wide <- reshape(
  summary,
  idvar = c("dataset", "target", "target_name"),
  timevar = "method",
  direction = "wide"
)
wide$clr_change_percent <- 100 * (
  wide$deviance.target_excluded_clr - wide$deviance.log1p
) / wide$deviance.log1p

# Reuse completed production benchmark scores for the expensive PLN comparison.
source(file.path(base_dir, "load_source.R"))
benchmark_env <- new.env()
load(file.path(base_dir, "bmr/hiv_lozupone_family.RData"), envir = benchmark_env)
benchmark_scores <- data.table::as.data.table(
  benchmark_env$bmr$score(poisson_measure)
)
benchmark_means <- benchmark_scores[
  learner_id %in% c("regr.featureless", "regr.pln"),
  .(deviance = mean(regr.poisson_deviance)),
  by = .(task_id, learner_id)
]
benchmark_wide <- reshape(
  as.data.frame(benchmark_means),
  idvar = "task_id", timevar = "learner_id", direction = "wide"
)
benchmark_wide$target_name <- sub("^Taxa", "", benchmark_wide$task_id)
names(benchmark_wide)[names(benchmark_wide) == "deviance.regr.featureless"] <- "deviance.featureless"
names(benchmark_wide)[names(benchmark_wide) == "deviance.regr.pln"] <- "deviance.pln"
wide <- merge(
  wide,
  benchmark_wide[, c("target_name", "deviance.featureless", "deviance.pln")],
  by = "target_name",
  all.x = TRUE,
  sort = FALSE
)

output_path <- file.path(base_dir, "out/pilot_clr_glmnet.csv")
write.csv(wide, output_path, row.names = FALSE)
print(wide, row.names = FALSE)
cat("Mean CLR change:", mean(wide$clr_change_percent), "%\n")
cat("Wrote", output_path, "\n")
