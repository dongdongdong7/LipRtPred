# This script is stored function used to build models.
# Barry Song
# 250224

# brnn: Bayesian Regularized Neural Networks
# rf: Random Forest
# xgboost: eXtreme Gradient Boosting

#' @rdname build_model
#' @title Build RT prediction model
#' @description
#' Build a machine learning model.
#'
#' @param trainingDf A training data frame.
#' @param k Either the number of folds or number of resampling iterations.
#' @param percentage For leave-group out cross-validation: the training percentage.
#' @param seed Seed number.
#' @param search Either "grid" or "random", describing how the tuning parameter grid is determined.
#' @param grid A data frame with possible tuning values.
#' @param metric A string that specifies what summary metric will be used to select the optimal model.
#' @param thread Number of threads in parallel.
#' @details
#' These functions use the caret packages.
#' For more information about parameters, please refer to \code{\link[caret]{train}} and \code{\link[caret]{trainControl}} .
#'
#' @return A model list.
#' @export
#'
#' @examples
#' model_rf <- build_rf(trainingDf = trainingDf)
build_rf <- function(trainingDf, k = 10, percentage = 0.8, seed = 1,
                     search = "random", grid = NULL, metric = "Rsquared",
                     thread = 1){
  control <- caret::trainControl(method = "cv",
                                 number = k,
                                 p = percentage,
                                 search = search,
                                 verboseIter = TRUE,
                                 allowParallel = TRUE)
  message("Building Random Forest model...")
  cl <- snow::makeCluster(thread)
  doSNOW::registerDoSNOW(cl)
  set.seed(seed)
  x <- trainingDf[, !colnames(trainingDf) %in% c("id", "smiles")]
  if(search == "grid"){
    if(is.null(grid)){
      grid <- base::expand.grid(mtry = c(1, 50, 100, 200, 500, 1000, 2000))
    }else grid <- grid
    model_rf <- caret::train(rt ~ .,
                             data = x,
                             method = "rf",
                             metric = metric,
                             trControl = control,
                             tuneGrid = grid)
  }else if(search == "random"){
    model_rf <- caret::train(rt ~ .,
                             data = x,
                             method = "rf",
                             metric = metric,
                             trControl = control)
  }
  snow::stopCluster(cl)
  gc()
  return(model_rf)
}
#' @rdname build_model
#' @export
#'
#' @examples
#' model_xgb <- build_xgb(trainingDf = trainingDf)
build_xgb <- function(trainingDf, k = 10, percentage = 0.8, seed = 1,
                      search = "random", grid = NULL, metric = "Rsquared",
                      thread = 1){
  cv.ctrl <- caret::trainControl(method = "cv", number = k, p = percentage,
                                 search = search, verboseIter = TRUE,
                                 allowParallel = TRUE)
  message("Building eXtreme Gradient Boosting model...")
  cl <- snow::makeCluster(thread)
  doSNOW::registerDoSNOW(cl)
  x <- trainingDf[, !colnames(trainingDf) %in% c("id", "smiles")]
  if(search == "grid" & is.null(grid)){
    if(is.null(gird)){
      # Use xgb.grid from Retip
      grid <- base::expand.grid(
        nrounds = c(300, 400, 500, 600, 700, 800, 1000),
        max_depth = c(2, 3, 4, 5),
        eta = c(0.01, 0.02),
        gamma = c(1),
        colsample_bytree = c(0.5),
        subsample = c(0.5),
        min_child_weight = c(10)
      )
    }
    else grid <- grid
    set.seed(seed)
    model_xgb <- caret::train(rt ~ .,
                              data = x,
                              method = "xgbTree",
                              metric = metric,
                              trControl = cv.ctrl,
                              tuneGrid = grid)
  }
  else if(search == "random"){
    set.seed(seed)
    model_xgb <- caret::train(rt ~ .,
                              data = x,
                              method = "xgbTree",
                              metric = metric,
                              trControl = cv.ctrl)
  }
  snow::stopCluster(cl)
  gc()
  return(model_xgb)
}
#' @rdname build_model
#' @export
#'
#' @examples
#' model_brnn <- build_brnn(trainingDf = trainingDf)
build_brnn <- function(trainingDf, k = 10, percentage = 0.8, seed = 1,
                       search = "random", grid = NULL, metric = "Rsquared",
                       thread = 1){
  # setting initial weight of neural network
  # seeds <- base::vector(mode = "list", length = nrow(trainingDf) + 1)
  # seeds <- base::lapply(seeds, function(x) 1:20)
  cv.ctrl <- caret::trainControl(method = "cv", number = k, p = percentage,
                                 search = search,
                                 verboseIter = TRUE,
                                 allowParallel = TRUE)
  message("Building Bayesian Regularized Neural Networks...")
  x <- trainingDf[, !colnames(trainingDf) %in% c("id", "smiles")]
  cl <- snow::makeCluster(thread)
  doSNOW::registerDoSNOW(cl)
  if(search == "grid"){
    if(is.null(grid)) tune.grid <- base::expand.grid(neurons = c(1, 2, 3, 4, 5))
    else tune.grid <- grid
    set.seed(seed)
    model_brnn <- caret::train(rt ~ ., data = x,
                               method = "brnn",
                               trControl = cv.ctrl,
                               metric = metric,
                               tuneGrid = tune.grid)
  }else if(search == "random"){
    set.seed(seed)
    model_brnn <- caret::train(rt ~ ., data = x,
                               method = "brnn",
                               trControl = cv.ctrl,
                               metric = metric)
  }
  snow::stopCluster(cl)
  gc()
  return(model_brnn)
}
#' @rdname build_model
#' @export
#'
#' @examples
#' model_glmnet <- build_glmnet(trainingDf = trainingDf)
build_glmnet <- function(trainingDf, k = 10, percentage = 0.8, seed = 1,
                         search = "random", grid = NULL, metric = "Rsquared",
                         thread = 1){
  cv.ctrl <- caret::trainControl(method = "cv", number = k, p = percentage,
                                 search = search,
                                 verboseIter = TRUE,
                                 allowParallel = TRUE)
  message("Building Elastic Net...")
  x <- trainingDf[, !colnames(trainingDf) %in% c("id", "smiles")]
  cl <- snow::makeCluster(thread)
  doSNOW::registerDoSNOW(cl)
  if(search == "grid"){
    if(is.null(grid)){
      tune.grid <- base::expand.grid(
        alpha = seq(0, 1, by = 0.1),
        lambda = 10^seq(-3, 3, length = 100)
      )
    }else tune.grid <- grid
    set.seed(seed)
    model_glmnet <- caret::train(rt ~ ., data = x,
                                 method = "glmnet",
                                 trControl = cv.ctrl,
                                 metric = metric,
                                 tuneGrid = tune.grid)
  }else if(search == "random"){
    set.seed(seed)
    model_glmnet <- caret::train(rt ~ ., data = x,
                                 method = "glmnet",
                                 trControl = cv.ctrl,
                                 metric = metric)
  }
  snow::stopCluster(cl)
  gc()
  return(model_glmnet)
}
