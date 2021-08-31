#' Compute metrics for model evaluation
#'
#' @param y_true A numeric vector with ground truth values.
#' @param y_pred A numeric vector with predicted values.
#'
#' @return A \code{data.frame} with computed metrics
compute_metrics <- function(y_true, y_pred){
    metrics <- tibble(
        corr = cor(y_true, y_pred),
        r2 = r2(y_true, y_pred),
        rse = rse(y_true, y_pred),
        sse = sse(y_true, y_pred),
        mse = mse(y_true, y_pred),
        rmse = rmse(y_true, y_pred),
        mae = mae(y_true, y_pred),
        mape = mape(y_true, y_pred)
    )
    return(metrics)
}


#' Sum of Squared Errors
#'
#' @param y_true A numeric vector with ground truth values.
#' @param y_pred A numeric vector with predicted values.
sse <- function(y_true, y_pred) {
    return(sum((y_true - y_pred) ** 2))
}


#' Sum of Squared Errors
#'
#' @param y_true A numeric vector with ground truth values.
#' @param y_pred A numeric vector with predicted values.
sse <- function(y_true, y_pred) {
    return(sum((y_true - y_pred) ** 2))
}


#' Relative Squared Error
#'
#' @param y_true A numeric vector with ground truth values.
#' @param y_pred A numeric vector with predicted values.
rse <- function(y_true, y_pred) {
    return(sse(y_true, y_pred) / sse(y_true, mean(y_true)))
}


#' \eqn{R^2} (coefficient of determination)
#'
#' @param y_true A numeric vector with ground truth values.
#' @param y_pred A numeric vector with predicted values.
r2 <- function(y_true, y_pred) {
    1 - rse(y_true, y_pred)
}


#' Mean Squared Error
#'
#' @param y_true A numeric vector with ground truth values.
#' @param y_pred A numeric vector with predicted values.
mse <- function(y_true, y_pred) {
    return(mean((y_true - y_pred) ** 2))
}


#' Root Mean Squared Error
#'
#' @param y_true A numeric vector with ground truth values.
#' @param y_pred A numeric vector with predicted values.
rmse <- function(y_true, y_pred) {
    return(sqrt(mse(y_true, y_pred)))
}


#' Absolute Error
#'
#' @param y_true A numeric vector with ground truth values.
#' @param y_pred A numeric vector with predicted values.
ae <- function(y_true, y_pred) {
    return(abs(y_true - y_pred))
}


#' Mean Absolute Error
#'
#' @param y_true A numeric vector with ground truth values.
#' @param y_pred A numeric vector with predicted values.
mae <- function(y_true, y_pred) {
    return(mean(ae(y_true, y_pred)))
}


#' Mean Absolute Percent Error
#'
#' @param y_true A numeric vector with ground truth values.
#' @param y_pred A numeric vector with predicted values.
mape <- function(y_true, y_pred) {
    return(mean(ae(y_true, y_pred) / abs(y_true)))
}


