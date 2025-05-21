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
  res <- stats::predict(model, x) # predcit return a numeric
  if(length(res) == nrow(x)){
    prd <- res
  }else{
    prd <- rep(NA, nrow(y))
    prd[as.integer(names(res))] <- res
  }

  prdDf <- data.frame(prd)
  colnames(prdDf) <- "rt_pred"
  resDf <- dplyr::as_tibble(cbind(y, prdDf)) %>%
    dplyr::filter(!is.na(rt_pred))
  return(resDf)
}
#' @rdname evaluate_model
#' @param digits integer indicating the number of decimal places (round) or significant digits (signif) to be used
#' @export
#' @examples
#' evaluate_model(testingDf = testingDf, model = model_rf)
evaluate_model <- function(testingDf, model, digits = 3){
  if(!any(colnames(testingDf) == "rt")){
    stop("testingDf do not have rt column!")
  }
  predDf <- predictRt(testingDf = testingDf, model = model)
  x <- testingDf[, !colnames(testingDf) %in% c("id", "smiles", "rt")]
  rt_max <- ceiling(max(max(predDf$rt), max(predDf$rt_pred))) + 1
  rmse <- round(sqrt(mean((predDf$rt - predDf$rt_pred)^2)), digits)
  mae <- round(mean(abs(predDf$rt - predDf$rt_pred)), digits)
  SSE <- sum((predDf$rt_pred - predDf$rt)^2)
  SST <- sum((mean(predDf$rt) - predDf$rt)^2)
  r_squared <- round(1- (SSE / SST), digits)
  predDf$difference <- round(abs(predDf$rt - predDf$rt_pred), digits = digits)
  max_point <- predDf %>%
    dplyr::slice_max(order_by = difference, n = 1)
  max_difference <- max(predDf$difference)
  n <- nrow(x)
  p <- ncol(x)
  if(n > p){
    r_squared_adjust <- round(1 - (1 - r_squared) * ((n - 1) / (n - p - 1)), 2)
    if(r_squared_adjust < 0) r_squared_adjust <- 0
  }
  else r_squared_adjust <- NA
  stats_table <- data.frame(Stats = c("MAE", "RMSE", "R2", "R2_adjust"), Values = c(mae, rmse, r_squared, r_squared_adjust))
  stats_table <- stats_table %>%
    dplyr::filter(!is.na(Values))
  stats_table <- gridExtra::tableGrob(stats_table, rows = NULL,
                                      cols = NULL, theme = gridExtra::ttheme_minimal(core = list(fg_params = list(hjust = 0,
                                                                                                                  x = 0.1)), base_colour = "#384049"))
  stats_table <- gtable::gtable_add_grob(stats_table, grobs = grid::rectGrob(gp = grid::gpar(fill = NA,
                                                                                             lwd = 2)), t = 1, b = nrow(stats_table), l = 1, r = ncol(stats_table))
  p1 <- ggplot2::ggplot(predDf, ggplot2::aes(rt, rt_pred)) +
    ggplot2::geom_point(ggplot2::aes(rt, rt_pred), colour = "#5E676F") +
    ggplot2::labs(title = paste0("Predicted vs Real - ", model$method), x="Observed RT", y = "Predicted RT") +
    ggplot2::xlim(0, rt_max) + ggplot2::ylim(0, rt_max) +
    ggplot2::annotation_custom(stats_table, xmin = rt_max * 0.125, xmax = rt_max * 0.250, ymin = rt_max * 0.750, ymax = rt_max  * 0.875) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(color = "#384049", face = "bold", hjust = 0.5), axis.line = ggplot2::element_line(colour = "#384049"), axis.text = ggplot2::element_text(colour = "#384049"), axis.title = ggplot2::element_text(colour = "#384049")) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "#D02937")
  p2 <- ggplot2::ggplot(predDf, ggplot2::aes(x = id, y = difference)) +
    ggplot2::geom_point(ggplot2::aes(x = id, y = difference), shape = 21, fill = "white", color = "black", alpha = 0.5) +
    ggrepel::geom_text_repel(mapping = ggplot2::aes(x = id,
                                                    y = difference,
                                                    label = paste0(id, "(", round(difference, digits), ")")),
                             data = max_point,
                             segment.color = "black",
                             min.segment.length = 0.5,
                             box.padding = 0.8) +
    ggplot2::labs(title = "", x = "", y = "Difference") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  return(list(predDf = predDf, MAE = mae, RMSE = rmse, R2 = r_squared, R2_adjust = r_squared_adjust, p1 = p1, Error_max = max_difference, p2 = p2))
}
