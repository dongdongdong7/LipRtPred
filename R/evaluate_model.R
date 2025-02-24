# Evaluate models using testing data.
# Barry Song
# 250224

#' @rdname evaluate_model
#' @title Evaluate model
#' @description
#' Evaluate model.
#'
#' @param testingDf A testing data frame.
#' @param model Model list.
#'
#' @export
#'
#' @examples
#' predictRt(testingDf = testingDf, model = model_rf)
predictRt <- function(testingDf, model){
  x <- testingDf[, !colnames(testingDf) %in% c("id", "smiles", "rt")]
  y <- testingDf[, colnames(testingDf) %in% c("id", "smiles", "rt")]
  prd <- data.frame(stats::predict(model, x))
  colnames(prd) <- "rt_pred"
  return(dplyr::as_tibble(cbind(y, prd)))
}
#' @rdname evaluate_model
#' @export
#' @examples
#' evaluate_model(testingDf = testingDf, model = model_rf)
evaluate_model <- function(testingDf, model){
  if(!any(colnames(testingDf) == "rt")){
    stop("testingDf do not have rt column!")
  }
  predDf <- predictRt(testingDf = testingDf, model = model)
  x <- testingDf[, !colnames(testingDf) %in% c("id", "smiles", "rt")]
  rt_max <- ceiling(max(max(predDf$rt), max(predDf$rt_pred))) + 1
  rmse <- round(sqrt(mean((predDf$rt - predDf$rt_pred)^2)), 2)
  mae <- round(mean(abs(predDf$rt - predDf$rt_pred)), 2)
  SSE <- sum((predDf$rt_pred - predDf$rt)^2)
  SST <- sum((mean(predDf$rt) - predDf$rt)^2)
  r_squared <- round(1- (SSE / SST), 2)
  n <- nrow(x)
  p <- ncol(x)
  if(n > p) r_squared_adjust <- round(1 - (1 - r_squared) * ((n - 1) / (n - p - 1)), 2)
  else r_squared_adjust <- NA
  stats_table <- data.frame(Stats = c("MAE", "RMSE", "R2", "R2_adjust"), Values = c(mae, rmse, r_squared, r_squared_adjust))
  stats_table <- stats_table %>%
    dplyr::filter(!is.na(Values))
  stats_table <- gridExtra::tableGrob(stats_table, rows = NULL,
                                      cols = NULL, theme = gridExtra::ttheme_minimal(core = list(fg_params = list(hjust = 0,
                                                                                                                  x = 0.1)), base_colour = "#384049"))
  stats_table <- gtable::gtable_add_grob(stats_table, grobs = grid::rectGrob(gp = grid::gpar(fill = NA,
                                                                                             lwd = 2)), t = 1, b = nrow(stats_table), l = 1, r = ncol(stats_table))
  p <- ggplot2::ggplot(predDf, ggplot2::aes(rt, rt_pred)) +
    ggplot2::geom_point(ggplot2::aes(rt, rt_pred), colour = "#5E676F") +
    ggplot2::labs(title = paste0("Predicted vs Real - ", model$method), x="Observed RT", y = "Predicted RT") +
    ggplot2::xlim(0, rt_max) + ggplot2::ylim(0, rt_max) +
    ggplot2::annotation_custom(stats_table, xmin = 1, xmax = 2, ymin = rt_max - 4, ymax = rt_max - 1) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(color = "#384049", face = "bold", hjust = 0.5), axis.line = ggplot2::element_line(colour = "#384049"), axis.text = ggplot2::element_text(colour = "#384049"), axis.title = ggplot2::element_text(colour = "#384049")) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "#D02937")
  return(list(MAE = mae, RMSE = rmse, R2 = r_squared, R2_adjust = r_squared_adjust, p = p))
}
