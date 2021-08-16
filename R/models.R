#' @import dplyr tibble
NULL


#' Fit (regularized) generalized linear model
#'
#' @param formula An object of class \code{formula} with a symbolic description
#' @param data A \code{data.frame} containing the variables in the model.
#' @param method A character string indicating the method to fit the model.
#' Possible values are \code{'glm'}, \code{'glmnet'} and \code{'cv.glmnet'}
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
    method = 'glm',
    family = gaussian,
    alpha = 1,
    ...
){
    result <- switch(
        method,
        'glm' = fit_glm(formula, data, family=family, ...),
        'glmnet' = fit_glmnet(formula, data, family=family, alpha=alpha, ...),
        'cv.glmnet' = fit_cvglmnet(formula, data, family=family, alpha=alpha, ...),
        'brms' = fit_brms(formula, data, family=family, ...),
        'xgb' = fit_xgb(formula, data, ...)
    )
    return(result)
}


#' Fit generalized linear model
#'
#' @param formula An object of class \code{formula} with a symbolic description
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
    formula, data,
    family = gaussian,
    alpha = 0.5,
    nlambda = 20,
    ...
){
    fit <- glmnetUtils::glmnet(
        formula,
        data = data,
        family = family,
        alpha = alpha,
        nlambda = nlambda,
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
#' @param data A \code{data.frame} containing the variables in the model.
#' @param family A description of the error distribution and link function to be used in the model.
#' See \code{\link[family]{stats}} for mode details.
#' @param alpha The elasticnet mixing parameter. See \code{\link[glmnet]{glmnet}} for details.
#' @param nlambda The number of \code{lambda} values.
#' @param nfolds The number of folds for CV.
#' @param ... Other parameters for the model fitting function.
#'
#' @return A list with two data frames: \code{gof} contains goodness of fit measures of the fit and
#' \code{coefs} contains the fitted coefficients.
#'
#' @export
fit_cvglmnet <- function(
    formula, data,
    family = gaussian,
    alpha = 0.5,
    nlambda = 20,
    nfolds = 5,
    ...
){
    fit <- glmnetUtils::cv.glmnet(
        formula,
        data = data,
        family = family,
        alpha = alpha,
        nlambda = nlambda,
        nfolds = nfolds,
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
    formula, data,
    family = gaussian,
    prior = brms::prior(normal(0,1)),
    ...
){
    # Silence annoying output
    sink('/dev/null', type='output')
    fit <- brms::brm(
        formula,
        data = data,
        family = family,
        prior = prior,
        silent = TRUE,
        refresh = 0,
        ...
    )
    sink()
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
#' @param data A \code{data.frame} containing the variables in the model.
#' @param ... Other parameters for the model fitting function.
#' @return A list with two data frames: \code{gof} contains goodness of fit measures of the fit and
#' \code{coefs} contains the fitted coefficients.
#'
#' @export
fit_xgb <- function(
    formula, data,
    params = list(
        max_depth=3,
        eta=0.01,
        objective='reg:squarederror'),
    nrounds = 1000,
    ...
){

    model_mat <- stats::model.matrix(formula, data=data)
    response <- data[[formula[[2]]]]
    fit <- xgboost::xgboost(
        data = model_mat,
        label = response,
        verbose = 0,
        params = params,
        nrounds = nrounds,
        ...
    )
    pred <- predict(fit, newdata=model_mat)
    resid <- sum((response - pred)**2)
    tot <- sum((response - mean(response))**2)
    gof <- tibble(
        rsq = 1 - resid / tot
    )
    coefs <- as_tibble(as.matrix(xgboost::xgb.importance(model=fit)))
    colnames(coefs) <- c('term', 'gain', 'cover', 'frequency')
    return(list(gof=gof, coefs=coefs))
}




