library(mlr3)
library(mlr3misc)
library(paradox)
library(PLNmodels)

pln_prepare_counts <- function(counts, offset_scheme, include_intercept = TRUE) {
  covariates <- if (include_intercept) {
    data.frame(Intercept = rep(1, nrow(counts)), row.names = rownames(counts))
  } else {
    data.frame(row.names = rownames(counts))
  }
  prepare_data(
    counts = counts,
    covariates = covariates,
    offset = offset_scheme %||% "TSS"
  )
}

pln_prepare_training_data <- function(task, offset_scheme, include_intercept = TRUE) {
  feature_names <- task$feature_names
  counts <- data.matrix(
    task$data(cols = c(task$target_names, feature_names)),
    rownames.force = TRUE
  )
  list(
    feature_names = feature_names,
    pln_data = pln_prepare_counts(counts, offset_scheme, include_intercept)
  )
}

pln_prepare_conditioning_inputs <- function(task, feature_names, offset_scheme) {
  feature_matrix <- data.matrix(task$data(cols = feature_names), rownames.force = TRUE)
  row_sums <- rowSums(feature_matrix)
  empty_samples <- which(row_sums == 0)
  non_empty_samples <- which(row_sums > 0)

  y_predict <- rep(NA_real_, nrow(feature_matrix))
  names(y_predict) <- rownames(feature_matrix)

  out <- list(
    empty_samples = empty_samples,
    non_empty_samples = non_empty_samples,
    y_predict = y_predict
  )

  if (length(non_empty_samples) > 0) {
    non_empty_counts <- feature_matrix[non_empty_samples, , drop = FALSE]
    test_pln_data <- pln_prepare_counts(non_empty_counts, offset_scheme)
    newdata <- data.frame(
      "(Intercept)" = rep(1, nrow(test_pln_data$Abundance)),
      check.names = FALSE,
      row.names = rownames(test_pln_data$Abundance)
    )
    out$test_pln_data <- test_pln_data
    out$newdata <- newdata
  }

  out
}

pln_finalize_predictions <- function(prep, y_predict, fallback_value) {
  if (length(prep$empty_samples) > 0) {
    y_predict[prep$empty_samples] <- 0
  }
  if (any(is.na(y_predict))) {
    y_predict[is.na(y_predict)] <- fallback_value
  }
  list(response = as.vector(y_predict))
}

pln_predict_with_conditioning <- function(task, feature_names, offset_scheme,
                                          model = NULL, fallback_value = NULL,
                                          prediction_is_valid = NULL,
                                          invalid_prediction_message = NULL) {
  prep <- pln_prepare_conditioning_inputs(task, feature_names, offset_scheme)
  y_predict <- prep$y_predict

  if (length(prep$non_empty_samples) > 0) {
    fallback_for_non_empty <- fallback_value %||% 0

    if (is.null(model)) {
      y_predict[prep$non_empty_samples] <- fallback_for_non_empty
    } else {
      valid_predictions <- predict_cond(
        model, prep$newdata, prep$test_pln_data$Abundance, type = "response"
      )[, 1]
      is_valid <- is.null(prediction_is_valid) || prediction_is_valid(valid_predictions)

      if (!is_valid) {
        if (!is.null(invalid_prediction_message)) warning(invalid_prediction_message)
        y_predict[prep$non_empty_samples] <- fallback_for_non_empty
      } else {
        y_predict[rownames(prep$test_pln_data$Abundance)] <- valid_predictions
      }
    }
  }

  pln_finalize_predictions(
    prep,
    y_predict,
    fallback_value %||% mean(y_predict, na.rm = TRUE)
  )
}

pln_network_predictions_are_valid <- function(valid_predictions) {
  !any(is.infinite(valid_predictions)) &&
    !any(is.nan(valid_predictions)) &&
    !any(valid_predictions < 0) &&
    !any(valid_predictions > 1e6) &&
    !(
      length(valid_predictions) > 1 &&
        max(valid_predictions, na.rm = TRUE) /
        min(valid_predictions[valid_predictions > 0], na.rm = TRUE) > 1e10
    )
}


pln_apply_shrinkage <- function(model, alpha) {
  sigma <- model$model_par$Sigma
  p     <- ncol(sigma)
  sigma[1, 2:p] <- alpha * sigma[1, 2:p]
  sigma[2:p, 1] <- alpha * sigma[2:p, 1]
  model$update(Sigma = sigma)
}

pln_poisson_deviance <- function(truth, response) {
  response  <- pmax(response, 1e-10)
  log_term  <- ifelse(truth == 0, 0, truth * log(truth / response))
  mean(2 * (log_term - (truth - response)))
}

LearnerRegrPLN <- R6::R6Class("LearnerRegrPLN",
  inherit = LearnerRegr,
  public = list(
    initialize = function() {
      ps <- ps(
        covariance  = p_fct(c("full", "diagonal", "spherical"), default = "full", tags = "train"),
        backend     = p_fct(c("nlopt", "torch"), default = "torch", tags = "train"),
        inner_folds = p_int(2L, 10L, default = 3L, tags = "train")
      )
      ps$values <- list(
        covariance = "full", backend = "torch", inner_folds = 3L
      )
      super$initialize(
        id = "regr.pln",
        param_set = ps,
        feature_types = c("integer", "numeric"),
        label = "Poisson Log-Normal Model",
        packages = "PLNmodels"
      )
    }
  ),
  private = list(
    .train = function(task) {
      pv           <- self$param_set$get_values(tags = "train")
      alpha_grid   <- seq(0.0, 1.0, by = 0.1)
      offset_scheme <- "none"
      n_inner      <- pv$inner_folds
      all_ids      <- task$row_ids
      n            <- length(all_ids)
      fold_ids     <- sample(rep_len(seq_len(n_inner), n))

      pln_fit <- function(t) {
        counts <- data.matrix(
          t$data(cols = c(t$target_names, t$feature_names)),
          rownames.force = TRUE
        )
        list(
          model         = PLN(
            Abundance ~ 1,
            data    = pln_prepare_counts(counts, offset_scheme),
            control = PLN_param(
              covariance = pv$covariance, trace = 0L, backend = pv$backend,
              config_optim = if (pv$backend == "torch") list(algorithm = "ADAM", lr = 0.01) else list()
            )
          ),
          feature_names = t$feature_names
        )
      }

      # Inner CV: fit PLN once per fold, sweep alpha cheaply
      alpha_pd <- numeric(length(alpha_grid))
      names(alpha_pd) <- as.character(alpha_grid)

      for (k in seq_len(n_inner)) {
        train_task <- task$clone(); train_task$filter(all_ids[fold_ids != k])
        val_task   <- task$clone(); val_task$filter(all_ids[fold_ids == k])
        truth      <- val_task$data(cols = task$target_names)[[1]]

        fit        <- pln_fit(train_task)
        sigma_orig <- fit$model$model_par$Sigma

        for (j in seq_along(alpha_grid)) {
          if (alpha_grid[j] < 1.0) pln_apply_shrinkage(fit$model, alpha_grid[j])
          pred <- pln_predict_with_conditioning(
            task          = val_task,
            feature_names = fit$feature_names,
            offset_scheme = offset_scheme,
            model         = fit$model
          )$response
          alpha_pd[j] <- alpha_pd[j] + pln_poisson_deviance(truth, pred) / n_inner
          if (alpha_grid[j] < 1.0) fit$model$update(Sigma = sigma_orig)
        }
      }

      best_alpha <- alpha_grid[which.min(alpha_pd)]

      # Final model on full training data with best alpha
      final <- pln_fit(task)
      if (best_alpha < 1.0) pln_apply_shrinkage(final$model, best_alpha)

      self$model <- list(
        pln_model     = final$model,
        feature_names = final$feature_names,
        best_alpha    = best_alpha,
        alpha_pd      = alpha_pd
      )
      invisible(self$model)
    },
    .predict = function(task) {
      pln_predict_with_conditioning(
        task          = task,
        feature_names = self$model$feature_names,
        offset_scheme = "none",
        model         = self$model$pln_model
      )
    }
  )
)

LearnerRegrPLNnetwork <- R6::R6Class("LearnerRegrPLNnetwork",
  inherit = LearnerRegr,
  public = list(
    initialize = function() {
      ps <- ps(
        rho = p_dbl(0.001, default = 0.1, tags = "train"),
        max_condition_number = p_dbl(1, 1e15, default = 1e12, tags = "train"),
        covariance = p_fct(c("full", "diagonal", "spherical"), default = "full", tags = "train"),
        trace = p_int(0, 2, default = 0, tags = "train"),
        backend = p_fct(c("nlopt", "torch"), default = "torch", tags = "train"),
        offset_scheme = p_fct(c("TSS", "CSS", "RLE", "GMPR", "none"), default = "none", tags = "train")
      )
      ps$values <- list(
        rho = 0.1,
        max_condition_number = 1e12,
        covariance = "full",
        trace = 0,
        backend = "torch",
        offset_scheme = "none"
      )
      super$initialize(
        id = "regr.pln_network",
        param_set = ps,
        feature_types = c("integer", "numeric"),
        label = "Poisson Log-Normal Network Model with Fallback",
        packages = "PLNmodels"
      )
    }
  ),
  private = list(
    .train = function(task) {
      pv <- self$param_set$get_values(tags = "train")
      prepared <- pln_prepare_training_data(task, pv$offset_scheme, include_intercept = FALSE)
      pln_data <- prepared$pln_data
      target_mean <- mean(task$data(cols = task$target_names)[[1]], na.rm = TRUE)

      inception <- PLN(
        Abundance ~ 1,
        data = pln_data,
        control = PLN_param(backend = pv$backend, covariance = pv$covariance, trace = pv$trace)
      )
      net_fit <- PLNnetwork(
        Abundance ~ 1,
        data = pln_data,
        penalties = pv$rho,
        control = PLNnetwork_param(
          backend = pv$backend,
          inception_cov = pv$covariance,
          trace = pv$trace,
          inception = inception
        )
      )

      best_model <- getBestModel(net_fit)
      network_condition_number <- tryCatch(
        kappa(sigma(best_model), exact = FALSE),
        error = function(...) Inf
      )
      if (!is.finite(network_condition_number) || network_condition_number > pv$max_condition_number) {
        warning("PLN network covariance is ill-conditioned, using featureless fallback")
        best_model <- NULL
      }

      self$model <- list(
        network_model = best_model,
        feature_names = prepared$feature_names,
        target_mean = target_mean,
        selected_penalty = if (is.null(best_model)) NA_real_ else best_model$penalty,
        network_density = if (is.null(best_model)) NA_real_ else best_model$density,
        network_condition_number = network_condition_number,
        pv = pv
      )
      invisible(self$model)
    },
    .predict = function(task) {
      featureless_prediction <- self$model$target_mean
      y_predict <- pln_predict_with_conditioning(
        task = task,
        feature_names = self$model$feature_names,
        offset_scheme = self$model$pv$offset_scheme,
        model = self$model$network_model,
        fallback_value = featureless_prediction,
        prediction_is_valid = pln_network_predictions_are_valid,
        invalid_prediction_message = "PLN network produced extreme predictions, using featureless fallback"
      )$response

      extreme_values <- is.infinite(y_predict) | is.nan(y_predict) | y_predict < 0 | y_predict > 1e6
      if (any(extreme_values)) y_predict[extreme_values] <- featureless_prediction

      list(response = as.vector(y_predict))
    }
  )
)

LearnerRegrPLNPCA <- R6::R6Class("LearnerRegrPLNPCA",
  inherit = LearnerRegr,
  public = list(
    initialize = function() {
      ps <- ps(
        trace = p_int(0, 2, default = 0, tags = "train"),
        backend = p_fct(c("nlopt", "torch"), default = "torch", tags = "train"),
        offset_scheme = p_fct(c("TSS", "CSS", "RLE", "GMPR", "none"), default = "none", tags = "train"),
        rank_max = p_int(1L, default = 5L, tags = "train"),
        criterion = p_fct(c("BIC", "ICL"), default = "BIC", tags = "train"),
        rho = p_dbl(0.001, default = 0.1, tags = "train")
      )
      ps$values <- list(
        trace = 0, backend = "torch", offset_scheme = "none",
        rank_max = 5L, criterion = "BIC", rho = 0.1
      )
      super$initialize(
        id = "regr.plnpca",
        param_set = ps,
        feature_types = c("integer", "numeric"),
        label = "Poisson Log-Normal PCA Model",
        packages = c("PLNmodels", "glassoFast")
      )
    }
  ),
  private = list(
    .train = function(task) {
      pv <- self$param_set$get_values(tags = "train")
      prepared <- pln_prepare_training_data(task, pv$offset_scheme)
      pln_data <- prepared$pln_data

      max_rank <- min(as.integer(pv$rank_max), nrow(pln_data$Abundance) - 1L, ncol(pln_data$Abundance) - 1L)
      if (max_rank < 1L) stop("PLNPCA requires at least two taxa and two non-empty training samples.")

      plnpca_family <- PLNPCA(
        Abundance ~ 1,
        data = pln_data,
        ranks = seq_len(max_rank),
        control = PLNPCA_param(
          backend = pv$backend, trace = pv$trace,
          config_optim = if (pv$backend == "torch") list(algorithm = "ADAM", lr = 0.01) else list()
        )
      )
      best_model <- getBestModel(plnpca_family, pv$criterion)
      pFix <- PLN(
        Abundance ~ 1,
        data = pln_data,
        control = PLN_param(
          covariance = "fixed",
          Omega = glassoFast::glassoFast(sigma(best_model), rho = pv$rho)$wi,
          backend = pv$backend,
          trace = pv$trace,
          config_optim = if (pv$backend == "torch") list(algorithm = "ADAM", lr = 0.01) else list()
        )
      )

      self$model <- list(
        pln_model = pFix,
        plnpca_family = plnpca_family,
        best_plnpca = best_model,
        feature_names = prepared$feature_names,
        pv = pv
      )
      invisible(self$model)
    },
    .predict = function(task) {
      pln_predict_with_conditioning(
        task = task,
        feature_names = self$model$feature_names,
        offset_scheme = self$model$pv$offset_scheme,
        model = self$model$pln_model
      )
    }
  )
)

LearnerRegrLasso <- R6::R6Class("LearnerRegrLasso",
  inherit = LearnerRegr,
  public = list(
    initialize = function() {
      super$initialize(
        id = "regr.lasso",
        feature_types = c("integer", "numeric"),
        label = "Poisson Lasso with Prediction Clipping",
        packages = "glmnet"
      )
    }
  ),
  private = list(
    .train = function(task) {
      x <- data.matrix(task$data(cols = task$feature_names))
      y <- task$data(cols = task$target_names)[[1]]
      target_mean <- mean(y, na.rm = TRUE)
      model <- tryCatch(
        glmnet::cv.glmnet(x, y, alpha = 1, family = "poisson", type.measure = "deviance"),
        error = function(e) NULL
      )
      self$model <- list(model = model, target_mean = target_mean, max_y = max(y, na.rm = TRUE))
      invisible(self$model)
    },
    .predict = function(task) {
      x <- data.matrix(task$data(cols = task$feature_names))
      fallback <- self$model$target_mean
      if (is.null(self$model$model)) {
        return(list(response = rep(fallback, nrow(x))))
      }
      response <- as.vector(
        predict(self$model$model, newx = x, s = "lambda.min", type = "response")
      )
      invalid <- is.infinite(response) | is.nan(response) | response < 0
      if (any(invalid)) response[invalid] <- fallback
      response <- pmin(response, self$model$max_y * 2)
      list(response = response)
    }
  )
)

MeasurePoissonDeviance <- R6::R6Class(
  "MeasurePoissonDeviance",
  inherit = mlr3::MeasureRegr,
  public = list(
    initialize = function() {
      super$initialize(
        id = "regr.poisson_deviance",
        range = c(0, Inf),
        minimize = TRUE,
        predict_type = "response",
        packages = character(0),
        properties = character(0),
        label = "Mean Poisson Deviance",
        man = "custom::poisson_deviance"
      )
    }
  ),
  private = list(
    .score = function(prediction, ...) {
      truth <- prediction$truth
      response <- prediction$response

      if (anyNA(response) || any(is.infinite(response))) {
        return(1e6)
      }

      response <- pmax(response, 1e-10)
      log_term <- ifelse(truth == 0, 0, truth * log(truth / response))
      mean(2 * (log_term - (truth - response)))
    }
  )
)


poisson_measure <- MeasurePoissonDeviance$new()
mlr3::mlr_measures$add("regr.poisson_deviance", MeasurePoissonDeviance)
