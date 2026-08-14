# CLR is computed from task features only. In each LOTO task, the target taxon
# is absent from these features and therefore cannot leak into the CLR center.
LearnerRegrCLRLasso <- R6::R6Class(
  "LearnerRegrCLRLasso",
  inherit = mlr3::LearnerRegr,
  public = list(
    initialize = function() {
      super$initialize(
        id = "regr.clr_lasso",
        feature_types = c("integer", "numeric"),
        label = "Target-excluded CLR Poisson Lasso",
        packages = "glmnet"
      )
    }
  ),
  private = list(
    .train = function(task) {
      x <- data.matrix(task$data(cols = task$feature_names))
      x <- x - rowMeans(x)
      y <- task$data(cols = task$target_names)[[1]]
      target_mean <- mean(y, na.rm = TRUE)
      model <- tryCatch(
        glmnet::cv.glmnet(
          x, y, alpha = 1, family = "poisson", type.measure = "deviance"
        ),
        error = function(e) NULL
      )
      self$model <- list(
        model = model,
        target_mean = target_mean,
        max_y = max(y, na.rm = TRUE),
        feature_names = task$feature_names
      )
      invisible(self$model)
    },
    .predict = function(task) {
      x <- data.matrix(task$data(cols = self$model$feature_names))
      x <- x - rowMeans(x)
      fallback <- self$model$target_mean
      if (is.null(self$model$model)) {
        return(list(response = rep(fallback, nrow(x))))
      }
      response <- as.vector(predict(
        self$model$model, newx = x, s = "lambda.min", type = "response"
      ))
      invalid <- !is.finite(response) | response < 0
      response[invalid] <- fallback
      response <- pmin(response, self$model$max_y * 2)
      list(response = response)
    }
  )
)
