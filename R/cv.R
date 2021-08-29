#' @import dplyr tibble
NULL


#' Cross-validation for GRN inference models
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' @param data A \code{data.frame} containing the variables in the model.
#' @param method A character string indicating the method to fit the model.
#' * \code{'glm'} - Generalized Liner Model with \code{\link[glm]{stats}}.
#' * \code{'glmnet'}, \code{'cv.glmnet'} - Regularized Generalized Liner Model with \code{\link[glmnet]{glmnet}}.
#' * \code{'brms'} - Bayesian Regression Models using \code{\link[brms-package]{brms}}.
#' * \code{'xgb'} - Gradient Boosting Regression using \code{\link[xgboost]{xgboost}}.
#' * \code{'bagging_ridge'} - Bagging Ridge Regression using scikit-learn via \link[xgboost]{reticulate}.
#' * \code{'bayesian_ridge'} - Bayesian Ridge Regression using scikit-learn via \link[xgboost]{reticulate}.
#' @param folds Number of cross-validation folds.
#' @param ... Other parameters for the model fitting function.
#'
#' @return A list with two data frames: \code{gof} contains goodness of fit measures of the fit and
#' \code{coefs} contains the fitted coefficients.
#'
#' @export
cv_model <- function(
    formula,
    data,
    method = c('glm', 'glmnet', 'cv.glmnet', 'brms', 'xgb', 'bagging_ridge', 'bayesian_ridge'),
    folds = 5,
    ...
){
    # Match args
    method <- match.arg(method)
    result <- switch(
        method,
        'glm' = cv_glm(formula, data, ...),
        'glmnet' = cv_glmnet(formula, data, ...),
        'cv.glmnet' = cv_cvglmnet(formula, data, ...),
        'brms' = cv_brms(formula, data, ...),
        'xgb' = cv_xgb(formula, data, ...),
        'bagging_ridge' = cv_bagging_ridge(formula, data, ...),
        'bayesian_ridge' = cv_bagging_ridge(formula, data, ...)
    )
    return(result)
}
