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
#' * \code{'xgb'} - Gradient Boosting Regression using \code{\link[xgboost]{xgboost}}.
#' * \code{'bagging_ridge'} - Bagging Ridge Regression using scikit-learn via \link[xgboost]{reticulate}.
#' * \code{'bayesian_ridge'} - Bayesian Ridge Regression using scikit-learn via \link[xgboost]{reticulate}.
#' @param k_folds Number of cross-validation folds.
#' @param strata Character vector with strata for stratified CV.
#' @param ... Other parameters for the model fitting function.
#'
#' @return A list with two data frames: \code{gof} contains goodness of fit measures of the fit and
#' \code{coefs} contains the fitted coefficients.
#'
#' @export
cv_model <- function(
    formula,
    data,
    method = c('glm', 'glmnet', 'cv.glmnet', 'xgb', 'bagging_ridge', 'bayesian_ridge'),
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
    y_pred <- suppressWarnings(predict(fit, test))
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
    which_max <- fit$index['1se', ]
    y_true <- test[[formula[[2]]]]
    y_pred <- as.numeric(predict(fit, newdata=test))
    metrics <- compute_metrics(y_true, y_pred)
    metrics$lambda <- fit$lambda.1se
    metrics$rsq <- fit$glmnet.fit$dev.ratio[which_max]
    metrics$alpha <- alpha
    return(metrics)
}


#' Score a gradient boosting regression model on a test set
#' using a range of regression metrics.
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' of the model to be fitted.
#' @param train A \code{data.frame} containing the training data.
#' @param test A \code{data.frame} containing the test data.
#' @param params A list with model parameters. For details, see \code{\link[xgb.train]{xgboost}}
#' @param ... Other parameters for the model fitting function.
#'
#' @return A \code{data.frame} containing different evaluation metrics.
#'
#' @export
score_xgb <- function(
    formula,
    train,
    test,
    params = list(
        max_depth=3,
        eta=0.01,
        objective='reg:squarederror'),
    nrounds = 1000,
    ...
){
    train_mat <- stats::model.matrix(formula, data=train)
    test_mat <- stats::model.matrix(formula, data=test)
    response <- train[[formula[[2]]]]
    fit <- xgboost::xgboost(
        data = train_mat,
        label = response,
        verbose = 0,
        params = params,
        nrounds = nrounds,
        ...
    )
    y_true <- test[[formula[[2]]]]
    y_pred <- predict(fit, newdata=test_mat)
    metrics <- compute_metrics(y_true, y_pred)
    metrics$rsq <- r2(response, predict(fit, newdata=train_mat))
    return(metrics)
}


#' Score a bagging ridge regression model on a test set
#' using a range of regression metrics.
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' of the model to be fitted.
#' @param train A \code{data.frame} containing the training data.
#' @param test A \code{data.frame} containing the test data.
#' @param alpha Positive float indicating the regularization strength.
#' @param solver Solver to use in the computational routines.
#' Options include ‘auto’, ‘svd’, ‘cholesky’, ‘lsqr’, ‘sparse_cg’, ‘sag’, ‘saga’.
#' @param bagging_number The number of ridge regression model in the bagging.
#' @param n_jobs The number of cores used to fit the model.
#' @param ... Other parameters for the model fitting function.
#'
#' @return A list with two data frames: \code{gof} contains goodness of fit measures of the fit and
#' \code{coefs} contains the fitted coefficients.
#'
#' @export
score_bagging_ridge <- function(
    formula,
    train,
    test,
    alpha = 1,
    solver = 'auto',
    bagging_number = 200L,
    n_jobs = -1,
    ...
){
    if (!require(reticulate, quietly = T)){
        stop('The reticulate package is required to use bagging ridge.')
    }
    np <- reticulate::import('numpy')
    pd <- reticulate::import('pandas')
    sklearn <- reticulate::import('sklearn')

    train_mat <- stats::model.matrix(formula, data=train)[,-1]
    test_mat <- stats::model.matrix(formula, data=test)[,-1]
    if (is.null(ncol(train_mat))){
        stop('The bagging ridge model requires at least two variables.')
    }
    response <- train[[formula[[2]]]]

    model <- sklearn$ensemble$BaggingRegressor(
        base_estimator = sklearn$linear_model$Ridge(
            alpha = alpha,
            solver = solver,
            random_state = as.integer(123),
            ...
        ),
        n_estimators = as.integer(bagging_number),
        bootstrap = TRUE,
        max_features = 0.8,
        n_jobs = as.integer(n_jobs),
        verbose = FALSE
    )
    model <- model$fit(train_mat, response)
    y_true <- test[[formula[[2]]]]
    y_pred <- model$predict(test_mat)
    metrics <- compute_metrics(y_true, y_pred)
    metrics$rsq <- r2(response, model$predict(train_mat))
    return(metrics)
}


#' Score a bagging ridge regression model on a test set
#' using a range of regression metrics.
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' of the model to be fitted.
#' @param train A \code{data.frame} containing the training data.
#' @param test A \code{data.frame} containing the test data.
#' @param ... Other parameters for the model fitting function.
#'
#' @return A list with two data frames: \code{gof} contains goodness of fit measures of the fit and
#' \code{coefs} contains the fitted coefficients.
#'
#' @export
score_bayesian_ridge <- function(
    formula,
    train,
    test,
    ...
){
    if (!require(reticulate, quietly = T)){
        stop('The reticulate package is required to use bayesian ridge.')
    }
    np <- reticulate::import('numpy')
    pd <- reticulate::import('pandas')
    sklearn <- reticulate::import('sklearn')

    model_mat <- stats::model.matrix(formula, data=data)[,-1]
    if (is.null(ncol(model_mat))){
        stop('The bayesian ridge model requires at least two variables.')
    }
    response <- train[[formula[[2]]]]

    model <- sklearn$linear_model$BayesianRidge(...)
    model <- model$fit(model_mat, response)
    model <- model$fit(train_mat, response)
    y_true <- test[[formula[[2]]]]
    y_pred <- model$predict(test_mat)
    metrics <- compute_metrics(y_true, y_pred)
    metrics$rsq <- model$score(model_mat, response)
    return(metrics)
}


