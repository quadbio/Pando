#' @import dplyr tibble
NULL


#' Fit (regularized) generalized linear model
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' of the model to be fitted.
#' @param data A \code{data.frame} containing the variables in the model.
#' @param method A character string indicating the method to fit the model.
#' * \code{'glm'} - Generalized Liner Model with \code{\link[glm]{stats}}.
#' * \code{'glmnet'}, \code{'cv.glmnet'} - Regularized Generalized Liner Model with \code{\link[glmnet]{glmnet}}.
#' * \code{'brms'} - Bayesian Regression Models using \code{\link[brms-package]{brms}}.
#' * \code{'xgb'} - Gradient Boosting Regression using \code{\link[xgboost]{xgboost}}.
#' * \code{'bagging_ridge'} - Bagging Ridge Regression using scikit-learn via \link[reticulate]{reticulate}.
#' * \code{'bayesian_ridge'} - Bayesian Ridge Regression using scikit-learn via \link[reticulate]{reticulate}.
#' @param family A description of the error distribution and link function to be used in the model.
#' See \code{\link[family]{stats}} for mode details.
#' @param alpha The elasticnet mixing parameter. See \code{\link[glmnet]{glmnet}} for details.
#' @param ... Other parameters for the model fitting function.
#'
#' @return A list with two data frames: \code{gof} contains goodness of fit measures of the fit and
#' \code{coefs} contains the fitted coefficients.
#'
#' @export
fit_model <- function(
    formula,
    data,
    method = c('glm', 'glmnet', 'cv.glmnet', 'brms', 'xgb', 'bagging_ridge', 'bayesian_ridge'),
    family = gaussian,
    alpha = 1,
    ...
){
    # Match args
    method <- match.arg(method)
    result <- switch(
        method,
        'glm' = fit_glm(formula, data, family=family, ...),
        'glmnet' = fit_glmnet(formula, data, family=family, alpha=alpha, ...),
        'cv.glmnet' = fit_cvglmnet(formula, data, family=family, alpha=alpha, ...),
        'brms' = fit_brms(formula, data, family=family, ...),
        'xgb' = fit_xgb(formula, data, ...),
        'bagging_ridge' = fit_bagging_ridge(formula, data, alpha=alpha, ...),
        'bayesian_ridge' = fit_bayesian_ridge(formula, data, ...)
    )
    return(result)
}


#' Fit generalized linear model
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' of the model to be fitted.
#' @param data A \code{data.frame} containing the variables in the model.
#' @param family A description of the error distribution and link function to be used in the model.
#' See \code{\link[family]{stats}} for mode details.
#' @param ... Other parameters for the model fitting function.
#'
#' @return A list with two data frames: \code{gof} contains goodness of fit measures of the fit and
#' \code{coefs} contains the fitted coefficients.
#'
#' @export
fit_glm <- function(formula, data, family=gaussian, ...){
    fit <- suppressWarnings(glm(formula, data=data, family=family, ...))
    s <- summary(fit)
    gof <- tibble(
        rsq = with(s, 1 - deviance/null.deviance)
    )
    coefs <- as_tibble(s$coefficients, rownames='term')
    colnames(coefs) <- c('term', 'estimate', 'std_err', 'statistic', 'pval')
    return(list(gof=gof, coefs=coefs))
}


#' Fit regularized generalized linear model
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' of the model to be fitted.
#' @param data A \code{data.frame} containing the variables in the model.
#' @param family A description of the error distribution and link function to be used in the model.
#' See \code{\link[family]{stats}} for mode details.
#' @param alpha The elasticnet mixing parameter. See \code{\link[glmnet]{glmnet}} for details.
#' @param nlambda The number of \code{lambda} values.
#' @param ... Other parameters for the model fitting function.
#'
#' @return A list with two data frames: \code{gof} contains goodness of fit measures of the fit and
#' \code{coefs} contains the fitted coefficients.
#'
#' @export
fit_glmnet <- function(
    formula,
    data,
    family = gaussian,
    alpha = 0.5,
    ...
){
    fit <- glmnetUtils::glmnet(
        formula,
        data = data,
        family = family,
        alpha = alpha,
        ...
    )
    class(fit) <- 'glmnet'
    which_max <- which(fit$dev.ratio > max(fit$dev.ratio) * 0.95)[1]
    lambda_choose <- fit$lambda[which_max]
    gof <- tibble(
        lambda = lambda_choose,
        rsq = fit$dev.ratio[which_max],
        alpha = alpha
    )
    coefs <- as_tibble(as.matrix(coef(fit, s=lambda_choose)), rownames='term')
    colnames(coefs) <- c('term', 'estimate')
    return(list(gof=gof, coefs=coefs))
}


#' Cross-validation for regularized generalized linear models
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' of the model to be fitted.
#' @param data A \code{data.frame} containing the variables in the model.
#' @param family A description of the error distribution and link function to be used in the model.
#' See \code{\link[family]{stats}} for mode details.
#' @param alpha The elasticnet mixing parameter. See \code{\link[glmnet]{glmnet}} for details.
#' @param ... Other parameters for the model fitting function.
#'
#' @return A list with two data frames: \code{gof} contains goodness of fit measures of the fit and
#' \code{coefs} contains the fitted coefficients.
#'
#' @export
fit_cvglmnet <- function(
    formula,
    data,
    family = gaussian,
    alpha = 0.5,
    ...
){
    fit <- glmnetUtils::cv.glmnet(
        formula,
        data = data,
        family = family,
        alpha = alpha,
        ...
    )
    class(fit) <- 'cv.glmnet'
    which_max <- fit$index['1se', ]
    gof <- tibble(
        lambda = fit$lambda.1se,
        rsq = fit$glmnet.fit$dev.ratio[which_max],
        alpha = alpha
    )
    coefs <- as_tibble(as.matrix(coef(fit)), rownames='term')
    colnames(coefs) <- c('term', 'estimate')
    return(list(gof=gof, coefs=coefs))
}


#' Fit a Bayesian regression model with brms and Stan
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' of the model to be fitted.
#' @param data A \code{data.frame} containing the variables in the model.
#' @param family A description of the error distribution and link function to be used in the model.
#' See \code{\link[family]{stats}} for mode details.
#' @param prior The prior distribution of the coefficients.
#' See \code{\link[set_prior]{brms}} for mode details.
#' The default (\code{prior(normal(0,1))}) results in ridge regularization.
#' @param ... Other parameters for the model fitting function.
#'
#' @return A list with two data frames: \code{gof} contains goodness of fit measures of the fit and
#' \code{coefs} contains the fitted coefficients.
#'
#' @export
fit_brms <- function(
    formula,
    data,
    family = gaussian,
    prior = brms::prior(normal(0,1)),
    ...
){
    if (!require(brms, quietly = T)){
        stop('The brms package is required to use brms bayesian regression models.')
    }
    fit <- suppressMessages(brms::brm(
        formula,
        data = data,
        family = family,
        prior = prior,
        silent = TRUE,
        refresh = 0,
        ...
    ))
    gof <- tibble(
        rsq = as.matrix(brms::bayes_R2(fit))[, 'Estimate']
    )
    coefs <- as_tibble(as.matrix(brms::fixef(fit, probs=c(0.05, 0.95))), rownames='term')
    colnames(coefs) <- c('term', 'estimate', 'est_error', 'q5', 'q95')
    coefs$pval <- bayestestR::p_map(fit)$p_MAP
    return(list(gof=gof, coefs=coefs))
}


#' Fit a gradient boosting regression model with XGBoost
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' of the model to be fitted.
#' @param data A \code{data.frame} containing the variables in the model.
#' @param params A list with model parameters. For details, see \code{\link[xgb.train]{xgboost}}
#' @param nrounds Max number of boosting iterations.
#' @param ... Other parameters for the model fitting function.
#'
#' @return A list with two data frames: \code{gof} contains goodness of fit measures of the fit and
#' \code{coefs} contains the fitted coefficients.
#'
#' @export
fit_xgb <- function(
    formula,
    data,
    params = list(
        max_depth=3,
        eta=0.01,
        objective='reg:squarederror'),
    nrounds = 1000,
    nthread = -1,
    ...
){
    if (!require(xgboost, quietly = T)){
        stop('The xgboost package is required to use xgb models.')
    }
    model_mat <- stats::model.matrix(formula, data=data)
    response <- data[[formula[[2]]]]
    fit <- xgboost::xgboost(
        data = model_mat,
        label = response,
        verbose = 0,
        params = params,
        nrounds = nrounds,
        nthread = nthread,
        ...
    )
    y_pred <- predict(fit, newdata=model_mat)
    gof <- tibble(
        rsq = r2(response, y_pred)
    )
    coefs <- as_tibble(as.data.frame(xgboost::xgb.importance(model=fit)))
    colnames(coefs) <- c('term', 'gain', 'cover', 'frequency')
    return(list(gof=gof, coefs=coefs))
}


#' Fit a bagging ridge regression model as implemented in scikit-learn (python)
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' of the model to be fitted.
#' @param data A \code{data.frame} containing the variables in the model.
#' @param alpha Positive float indicating the regularization strength.
#' @param solver Solver to use in the computational routines.
#' Options include ‘auto’, ‘svd’, ‘cholesky’, ‘lsqr’, ‘sparse_cg’, ‘sag’, ‘saga’.
#' @param bagging_number The number of ridge regression model in the bagging.
#' @param n_jobs The number of cores used to fit the model.
#' @param p_method The test used to calculate p-values. Options are 't' for \code{t.test}, and 'wilcox' for \code{wilcox.test}
#' @param ... Other parameters for the model fitting function.
#'
#' @return A list with two data frames: \code{gof} contains goodness of fit measures of the fit and
#' \code{coefs} contains the fitted coefficients.
#'
#' @export
fit_bagging_ridge <- function(
    formula,
    data,
    alpha = 1,
    solver = 'auto',
    bagging_number = 200L,
    n_jobs = 1,
    p_method = c('wilcox', 't'),
    ...
){
    p_method <- match.arg(p_method)
    if (!require(reticulate, quietly = T)){
        stop('The reticulate package is required to use bagging ridge models.')
    }
    np <- reticulate::import('numpy')
    pd <- reticulate::import('pandas')
    sklearn <- reticulate::import('sklearn')

    model_mat <- stats::model.matrix(formula, data=data)[,-1]
    if (is.null(ncol(model_mat))){
        stop('The bagging ridge model requires at least two variables.')
    }
    response <- data[[formula[[2]]]]

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
    model <- model$fit(model_mat, response)

    idx_features <- do.call(cbind, model$estimators_features_) + 1
    coefs_features <- sapply(model$estimators_, function(x) x$coef_)
    coefs <- t(sapply(1:bagging_number, function(i){
        coefs <- setNames(rep(NaN, ncol(model_mat)), colnames(model_mat))
        if (ncol(model_mat) > 2){
            coefs[idx_features[,i]] <- coefs_features[,i]
        }
        if (ncol(model_mat) == 2){
            coefs[idx_features[,i]] <- coefs_features[i]
        }
        return(coefs)
    }))
    p <- switch(
        p_method,
        't' = apply(coefs, 2, function(x) t.test(x[!is.nan(x)])$p.value),
        'wilcox' = apply(coefs, 2, function(x) wilcox.test(x[!is.nan(x)])$p.value)
    )
    coefs <- tibble(
        term = colnames(model_mat),
        estimate = colMeans(coefs, na.rm=TRUE),
        pval = p,
        neglog10p = -log10(ifelse(is.na(p), 1, p))
    )
    y_pred <- model_mat %*% matrix(coefs$estimate)
    gof <- tibble(
        rsq = r2(response, y_pred)
    )
    return(list(gof=gof, coefs=coefs))
}


#' Fit a bayesian ridge regression model as implemented in scikit-learn (python)
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' of the model to be fitted.
#' @param data A \code{data.frame} containing the variables in the model.
#' @param ... Other parameters for the model fitting function.
#'
#' @return A list with two data frames: \code{gof} contains goodness of fit measures of the fit and
#' \code{coefs} contains the fitted coefficients.
#'
#' @export
fit_bayesian_ridge <- function(
    formula,
    data,
    ...
){
    if (!require(reticulate, quietly = T)){
        stop('The reticulate package is required to use bayesian ridge models.')
    }
    np <- reticulate::import('numpy')
    pd <- reticulate::import('pandas')
    sklearn <- reticulate::import('sklearn')

    model_mat <- stats::model.matrix(formula, data=data)[,-1]
    if (is.null(ncol(model_mat))){
        stop('The bayesian ridge model requires at least two variables.')
    }
    response <- data[[formula[[2]]]]

    model <- sklearn$linear_model$BayesianRidge(...)
    model <- model$fit(model_mat, response)

    coefs <- model$coef_
    coef_var <- diag(model$sigma_)
    p <- pnorm(q=0, mean=abs(coefs), sd=sqrt(coef_var)) * 2
    coefs <- tibble(
        term = colnames(model_mat),
        estimate = coefs,
        est_variance = coef_var,
        pval = p,
        neglog10p = -log10(ifelse(is.na(p), 1, p))
    )
    gof <- tibble(
        rsq = model$score(model_mat, response)
    )
    return(list(gof=gof, coefs=coefs))
}

#' \eqn{R^2} (coefficient of determination)
#'
#' @param y_true A numeric vector with ground truth values.
#' @param y_pred A numeric vector with predicted values.
r2 <- function(y_true, y_pred) {
    1 - rse(y_true, y_pred)
}

#' Relative Squared Error
#'
#' @param y_true A numeric vector with ground truth values.
#' @param y_pred A numeric vector with predicted values.
rse <- function(y_true, y_pred) {
    return(sse(y_true, y_pred) / sse(y_true, mean(y_true)))
}

#' Sum of Squared Errors
#'
#' @param y_true A numeric vector with ground truth values.
#' @param y_pred A numeric vector with predicted values.
sse <- function(y_true, y_pred) {
    return(sum((y_true - y_pred) ** 2))
}



