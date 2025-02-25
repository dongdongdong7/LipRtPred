# Preprocess data by dividing it into test and training sets and removing redundant descriptors.
# Barry Song
# 250223

#' @rdname preprocess_data
#' @title Preprocess raw data
#' @description
#' Preprocess raw data as follows: (1) remove columns containing NA; (2) remove columns with variance close to 0.
#' The data were then divided into training and testing data.
#'
#' @param inputDf Input data frame with id, smiles and MD or FP.
#' @param freqCut The cutoff for the ratio of the most common value to the second most common value.
#' @param uniqueCut The cutoff for the percentage of distinct values out of the number of total samples.
#' @param thread Number of threads in parallel.
#' @param cutoff A numeric value for the pair-wise absolute correlation cutoff.
#' @details
#' This function uses the caret package, for more information on \code{\link[caret]{nearZeroVar}}.
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' data("cmpDf_demo", package = "LipRtPred")
#' descsDf <- GetCDK_MD(cmpDf = cmpDf_demo)
#' inputDf <- filterColumns(inputDf = descsDf)
filterColumns <- function(inputDf, freqCut = 95/5, uniqueCut = 10, thread = 1, cutoff = 0.9){
  marksColumns <- which(colnames(inputDf) %in% c("id", "smiles", "rt"))
  marksDf <- inputDf[, marksColumns]
  columnsDf <- inputDf[, -marksColumns]

  message("Remove NA columns...")
  columnsDf <- columnsDf[, !apply(columnsDf, 2, function(x) any(is.na(x)))]

  cl <- snow::makeCluster(thread)
  doSNOW::registerDoSNOW(cl)
  message("Remove columns with near zero variance values...")
  nzvColumns <- caret::nearZeroVar(columnsDf, saveMetrics = FALSE, freqCut = freqCut, uniqueCut = uniqueCut, allowParallel = TRUE)
  if(length(nzvColumns) != 0) columnsDf <- columnsDf[, -nzvColumns]
  snow::stopCluster(cl)
  gc()

  message("Remove columns with high correlation...")
  R <- stats::cor(columnsDf)
  corColumns <- caret::findCorrelation(R, cutoff = cutoff)
  if(length(corColumns) != 0) columnsDf <- columnsDf[, -corColumns]

  return(dplyr::as_tibble(cbind(marksDf, columnsDf)))
}
#' @rdname preprocess_data
#' @param percentage The percentage of data that goes to training.
#' @param seed The seed number.
#'
#' @return A list with training and testing sets.
#' @export
#'
#' @examples
#' tmp <- split_data(inputDf = inputDf, percentage = 0.8)
#' trainingDf <- tmp$training
#' testingDf <- tmp$testing
split_data <- function(inputDf, percentage = 0.8, seed = 1){
  set.seed(seed)
  idx_training <- as.integer(unlist(caret::createDataPartition(1:nrow(inputDf), p = percentage, list = TRUE)))
  return(list(training = inputDf[idx_training, ], testing = inputDf[-idx_training, ]))
}

