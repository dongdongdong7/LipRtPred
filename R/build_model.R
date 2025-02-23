# .build_rf(trainingDf = trainingDf)
.build_rf <- function(trainingDf, k = 10, search = "random", seed = 1, thread = 1,
                      metric = "Rsquared"){
  browser()
  control <- caret::trainControl(method = "cv",
                                 number = k,
                                 search = search,
                                 verboseIter = TRUE)
  message("Building Random Forest model...")
  set.seed(seed)
  x <- trainingDf[, !colnames(trainingDf) %in% c("id", "smiles")]
  cores <- parallel::makeCluster(thread)
  doParallel::registerDoParallel(cores = cores)
  if(thread == 1){
    model_rf <- caret::train(rt ~ .,
                             data = x,
                             method = "rf",
                             metric = metric,
                             trControl = control,
                             importance = TRUE,
                             allowParallel = FALSE)
  }else{
    model_rf <- caret::train(rt ~ .,
                             data = x,
                             method = "rf",
                             metric = metric,
                             trControl = control,
                             importance = TRUE,
                             allowParallel = TRUE)
  }
  parallel::stopCluster(cores)
  gc()
  return(model_rf)
}
