#' @import dplyr tibble
NULL


#' Cross-validation for GRN inference models
#'
#' @importFrom purrr map_dfr set_names
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' of the model to be fitted.
#' @param data A \code{data.frame} containing the variables in the model.
#' @param method A character string indicating the method to fit the model.
#' * \code{'glm'} - Generalized Liner Model with \code{\link[glm]{stats}}.
#' * \code{'glmnet'}, \code{'cv.glmnet'} - Regularized Generalized Liner Model with \code{\link[glmnet]{glmnet}}.
#' * \code{'brms'} - Bayesian Regression Models using \code{\link[brms-package]{brms}}.
#' * \code{'xgb'} - Gradient Boosting Regression using \code{\link[xgboost]{xgboost}}.
#' * \code{'bagging_ridge'} - Bagging Ridge Regression using scikit-learn via \link[xgboost]{reticulate}.
#' * \code{'bayesian_ridge'} - Bayesian Ridge Regression using scikit-learn via \link[xgboost]{reticulate}.
#' @param k_folds Number of cross-validation folds.
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
    k_folds = 5,
    strata = NULL,
    ...
){
    # Match args
    method <- match.arg(method)
    # Get scoring function
    score_func <- switch(
        method,
        'glm' = score_glm,
        'glmnet' = score_glmnet,
        'cv.glmnet' = score_cvglmnet,
        'brms' = score_brms,
        'xgb' = score_xgb,
        'bagging_ridge' = score_bagging_ridge,
        'bayesian_ridge' = score_bagging_ridge
    )
    # Loop through folds
    folds <- cv_folds(data, folds=k_folds, strata=strata)
    result <- map_dfr(set_names(folds, seq_along(folds)), function(idx){
        score_func(formula, train=data[-idx, ], test=data[idx, ], ...)
    }, .id='fold')
    return(result)
}


#' Score a generalized linear model on a test set using a range of regression metrics.
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' of the model to be fitted.
#' @param train A \code{data.frame} containing the training data.
#' @param test A \code{data.frame} containing the test data.
#' @param ... Other parameters for the model fitting function.
#'
#' @return A \code{data.frame} containing different evaluation metrics.
#'
#' @export
score_glm <- function(formula, train, test, ...){
    fit <- suppressWarnings(glm(formula, data=train, ...))
    s <- summary(fit)
    y_true <- test[[formula[[2]]]]
    y_pred <- predict(fit, test)
    metrics <- compute_metrics(y_true, y_pred)
    metrics$dsq <- with(s, 1 - deviance/null.deviance)
    return(metrics)
}


#' Score a regularized generalized linear model on a test set using a range of regression metrics.
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' of the model to be fitted.
#' @param train A \code{data.frame} containing the training data.
#' @param test A \code{data.frame} containing the test data.
#' @param alpha The elasticnet mixing parameter. See \code{\link[glmnet]{glmnet}} for details.
#' @param ... Other parameters for the model fitting function.
#'
#' @return A \code{data.frame} containing different evaluation metrics.
#'
#' @export
score_glmnet <- function(formula, train, test, alpha=0.5, ...){
    fit <- glmnetUtils::glmnet(
        formula,
        data = train,
        alpha = alpha,
        ...
    )
    which_max <- which(fit$dev.ratio > max(fit$dev.ratio) * 0.95)[1]
    lambda_choose <- fit$lambda[which_max]
    y_true <- test[[formula[[2]]]]
    y_pred <- predict(fit, newdata=test)[, which_max]
    metrics <- compute_metrics(y_true, y_pred)
    metrics$lambda <- lambda_choose
    metrics$dsq <- fit$dev.ratio[which_max]
    metrics$alpha <- alpha
    return(metrics)
}


#' Score a cross-validated regularized generalized linear model on a test set
#' using a range of regression metrics.
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' of the model to be fitted.
#' @param train A \code{data.frame} containing the training data.
#' @param test A \code{data.frame} containing the test data.
#' @param alpha The elasticnet mixing parameter. See \code{\link[glmnet]{glmnet}} for details.
#' @param ... Other parameters for the model fitting function.
#'
#' @return A \code{data.frame} containing different evaluation metrics.
#'
#' @export
score_cvglmnet <- function(formula, train, test, alpha=0.5, ...){
    fit <- glmnetUtils::cv.glmnet(
        formula,
        data = train,
        alpha = alpha,
        ...
    )
    class(fit) <- 'cv.glmnet'
    which_max <- fit$index['1se', ]
    y_true <- test[[formula[[2]]]]
    y_pred <- predict(fit, newdata=test)[, which_max]
    metrics <- compute_metrics(y_true, y_pred)
    metrics$lambda <- fit$lambda.1se,
    metrics$rsq <- fit$glmnet.fit$dev.ratio[which_max],
    metrics$alpha <- alpha
    return(metrics)
}







