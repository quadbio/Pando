library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Pando)
library(doParallel)
library(devtools)
library(Signac)
library(Seurat)
library(brms)
library(bayestestR)

rename <- dplyr::rename

data(motifs)
data(phastConsElements20Mammals.UCSC.hg38)
data(EnsDb.Hsapiens.v93.annot.UCSC.hg38)
data(motif2tf)

genes_use <- c('PAX6', 'POU5F1', 'DLX6', 'NFIA')
tfs_use <- c('PAX6', 'POU5F1', 'DLX6', 'NFIA')

motif2tf_use <- motif2tf %>%
    filter(origin=='JASPAR2020') %>%
    rename('mot'=motif)

registerDoParallel(2)

test_srt <- read_rds('../data/test_seurat.rds')
test_srt <- initiate_grn(
    test_srt,
    genes = genes_use,
    regions = phastConsElements20Mammals.UCSC.hg38
)
test_srt <- find_motifs(
    test_srt,
    motif_tfs = motif2tf_use,
    pfm = motifs[1:10],
    genome = BSgenome.Hsapiens.UCSC.hg38
)
test_srt <- infer_grn(test_srt, peak_to_gene_method = 'Signac', parallel=T)
test_srt <- infer_grn(test_srt, peak_to_gene_method = 'GREAT', parallel=T, method = 'cv.glmnet', nlambda=100, alpha=0.3)
test_srt <- infer_grn(test_srt, peak_to_gene_method = 'GREAT', method = 'brms', backend='cmdstanr', prior=prior(lasso()))

test_srt <- find_modules(test_srt, min_genes_per_module=0)

NetworkParams(test_srt)$method

coef(test_srt)
gof(test_srt)
modules <- NetworkModules(test_srt)

expect_true('tbl'%in%class(modules@meta))

class(modules@meta)


### Try random forest (xgboost) ####

np <- 100
x <- seq(np) + rnorm(np, sd=3)
z <- seq(np) + rnorm(np, sd=5)
y <- round(rnorm(np, sd=50, mean=1000) + x * 3 + z * 1.5)

tbl <- tibble(
    x = x,
    y = y,
    z = z
)

xy <- model_as_xy(data=as.data.frame(tbl), formula=y ~ x:z)


formula <- y ~ x:z
data <- tbl
model_mat <- stats::model.matrix(formula, data = data)
response <- data[[formula[[2]]]]


params <- list(
    max_depth = 3
)

object <- xgboost::xgboost(data = model_mat, label = response, nrounds=100, params=params)
pred <- predict(object, newdata=model_mat)
resid_sq <- (pred - response)**2

sstot <- sum((pred - mean(response))^2)
ssresid <- sum(resid_sq)
1 - ssresid / sstot

attr(object, 'formula') <- formula

summary(object)
xgboost::xgb.importance(model=object)















