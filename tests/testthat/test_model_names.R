context('Tests for placeholder variable names in model formulas')

library(testthat)
library(Pando)


# Names which are not valid R symbols: reserved words, parentheses, colons,
# leading digits (e.g. dm6 chromosomes '3R') and underscores
target <- 'if'
regions <- c('3R-100-200', '2L-5-50', 'chr1-1-9')
tfs <- c('Su(H)', 'lncRNA:CR1', '18w', 'Abd-B', 'ey_1')

set.seed(1)
n <- 60
gene_data <- matrix(rnorm(n*7), n, dimnames=list(NULL, c(target, tfs, 'other')))
peak_x <- matrix(runif(n*3), n, dimnames=list(NULL, regions))
gene_data[, target] <- 2*peak_x[, 1]*gene_data[, 'Su(H)'] + rnorm(n, sd=0.1)

# 'Su(H)' is shared by two regions, so glm orders some terms as 'tf:region'
peak_tfs <- list(
    `3R-100-200` = c('Su(H)', 'lncRNA:CR1'),
    `2L-5-50` = c('18w', 'Abd-B', 'ey_1', 'Su(H)')
)
pairs <- c(
    '3R-100-200 Su(H)', '3R-100-200 lncRNA:CR1', '2L-5-50 18w',
    '2L-5-50 Abd-B', '2L-5-50 ey_1', '2L-5-50 Su(H)'
)

fit_and_decode <- function(method, term, peak_tfs_use=peak_tfs){
    input <- Pando:::build_model_input(target, peak_tfs_use, gene_data, peak_x, term=term)
    fit <- Pando:::fit_model(input$formula, data=input$data, method=method)
    coefs <- Pando:::decode_model_terms(fit$coefs, input$var_names)
    format_coefs(cbind(target=target, coefs), term=term)
}

methods <- c('glm', 'glmnet', 'cv.glmnet', 'xgb')
method_pkgs <- c(glm='stats', glmnet='glmnet', cv.glmnet='glmnet', xgb='xgboost')


test_that('model_var_names gives unique syntactic names.', {
    var_names <- Pando:::model_var_names(target, regions, tfs)
    expect_equal(names(var_names), c(target, regions, tfs))
    expect_equal(unname(make.names(var_names)), unname(var_names))
    expect_false(any(duplicated(var_names)))
})

test_that('model_var_names fails on a gene named like a region.', {
    expect_error(Pando:::model_var_names('chr1-5-9', 'chr1-5-9', 'A'), 'chr1-5-9')
})

test_that('build_model_input makes a valid formula and model matrix.', {
    input <- Pando:::build_model_input(target, peak_tfs, gene_data, peak_x, term=':')
    expect_equal(input$nfeats, 6)
    expect_equal(all.vars(input$formula)[1], 'target')
    expect_true(all(all.vars(input$formula) %in% colnames(input$data)))
    expect_false(anyNA(colnames(input$data)))
    # All regions are in the model matrix, genes without TF role are not
    expect_equal(ncol(input$data), 1 + 5 + 3)
    expect_false('other' %in% names(input$var_names))
    # Data columns keep their values under the new names
    expect_equal(input$data$target, unname(gene_data[, target]))
    expect_equal(input$data[[input$var_names['3R-100-200']]], unname(peak_x[, '3R-100-200']))
})

test_that('decode_model_terms maps placeholders back to original names.', {
    var_names <- Pando:::model_var_names(target, regions, c(tfs, paste0('X', 1:5)))
    coefs <- data.frame(
        term = c('(Intercept)', 'region1:tf1', 'tf2:region1', 'region2', 'tf5', 'region2:tf10'),
        estimate = 1:6
    )
    dec <- Pando:::decode_model_terms(coefs, var_names)
    expect_equal(dec$tf, c(NA, 'Su(H)', 'lncRNA:CR1', NA, 'ey_1', 'X5'))
    expect_equal(dec$region, c(NA, '3R-100-200', '3R-100-200', '2L-5-50', NA, '2L-5-50'))
    expect_equal(dec$term[2:3], c('3R-100-200:Su(H)', 'lncRNA:CR1:3R-100-200'))
})

test_that('decode_model_terms keeps unknown placeholders in term.', {
    var_names <- Pando:::model_var_names(target, regions[1], tfs[1])
    dec <- Pando:::decode_model_terms(data.frame(term='region1:tf9'), var_names)
    expect_equal(dec$term, '3R-100-200:tf9')
    expect_true(is.na(dec$tf))
})

test_that('interaction models decode to the expected region-TF pairs.', {
    for (method in methods){
        skip_if_not_installed(method_pkgs[[method]])
        coefs <- fit_and_decode(method, ':')
        # xgb reports feature gain instead of an estimate
        coefs$value <- if ('estimate' %in% colnames(coefs)) coefs$estimate else coefs$gain
        coefs <- coefs[coefs$value != 0, ]
        expect_false(any(c('(Intercept)', 'Intercept') %in% coefs$term), info=method)
        expect_true(all(coefs$target == target), info=method)
        expect_false(anyNA(coefs$tf) || anyNA(coefs$region), info=method)
        expect_true(all(paste(coefs$region, coefs$tf) %in% pairs), info=method)
        if (method == 'glm'){
            expect_setequal(paste(coefs$region, coefs$tf), pairs)
        }
        # The true effect is recovered under its original names
        top <- coefs[which.max(abs(coefs$value)), ]
        expect_equal(c(top$region, top$tf), c('3R-100-200', 'Su(H)'), info=method)
    }
})

test_that('main effects and interactions are split correctly for "*".', {
    for (method in c('glm', 'xgb')){
        skip_if_not_installed(method_pkgs[[method]])
        coefs <- fit_and_decode(method, '*')
        inter <- grepl(':', gsub('lncRNA:CR1', '', coefs$term, fixed=TRUE))
        expect_false(anyNA(coefs$tf[inter]) || anyNA(coefs$region[inter]), info=method)
        expect_true(all(xor(is.na(coefs$tf[!inter]), is.na(coefs$region[!inter]))), info=method)
        expect_true(all(coefs$tf %in% c(tfs, NA)), info=method)
        expect_true(all(coefs$region %in% c(regions, NA)), info=method)
    }
})

test_that('additive models give rows with either TF or region.', {
    coefs <- fit_and_decode('glm', '+')
    expect_true(all(xor(is.na(coefs$tf), is.na(coefs$region))))
    expect_setequal(na.omit(coefs$tf), c('Su(H)', 'lncRNA:CR1', '18w', 'Abd-B', 'ey_1'))
    expect_setequal(na.omit(coefs$region), c('3R-100-200', '2L-5-50'))
})

test_that('a single region with a single TF works.', {
    coefs <- fit_and_decode('glm', ':', peak_tfs_use=list(`3R-100-200`='Su(H)'))
    expect_equal(nrow(coefs), 1)
    expect_equal(c(coefs$tf, coefs$region, coefs$term), c('Su(H)', '3R-100-200', '3R-100-200:Su(H)'))
})

test_that('format_coefs still parses terms without tf/region columns.', {
    coefs <- data.frame(
        target = 'Gata',
        term = c('(Intercept)', 'chr1_1_9:Abd_B'),
        estimate = c(0, 1)
    )
    res <- format_coefs(coefs)
    expect_equal(res$tf, 'Abd-B')
    expect_equal(res$region, 'chr1-1-9')
})
