library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(doParallel)
library(devtools)
library(Signac)
library(Seurat)
library(brms)
library(bayestestR)
library(Pando)

rename <- dplyr::rename

data(motifs)
data(phastConsElements20Mammals.UCSC.hg38)
data(SCREEN.ccRE.UCSC.hg38)
data(motif2tf)

genes_use <- c('PAX6', 'POU5F1', 'DLX6', 'NFIA')
tfs_use <- c('PAX6', 'POU5F1', 'DLX6', 'NFIA')

motif2tf_use <- motif2tf %>%
    filter(origin=='JASPAR2020') %>%
    rename('mot'=motif)

registerDoParallel(2)

test_srt <- read_rds('../data/test_seurat.rds')

test_mat <- test_srt[['RNA']]@data
rownames(test_mat)[rownames(test_mat)=='NFIA'] <- '2610316D01Rik'
test_srt2 <- CreateSeuratObject(test_mat, assay='RNA')
test_srt2[['peaks']] <- test_srt[['peaks']]
motif2tf_use$tf[motif2tf_use$tf=='NFIA'] <- '2610316D01Rik'

test_srt <- initiate_grn(
    test_srt2,
    regions = phastConsElements20Mammals.UCSC.hg38
)

test_srt <- find_motifs(
    test_srt,
    motif_tfs = motif2tf_use,
    pfm = motifs[1:10],
    genome = BSgenome.Hsapiens.UCSC.hg38
)

NetworkTFs(test_srt)

test_srt <- aggregate_assay(test_srt, group_name = 'seurat_clusters')

test_srt <- infer_grn(test_srt, genes=genes_use,
    peak_to_gene_method = 'GREAT', parallel=F,
    aggregate_peaks_col='seurat_clusters', aggregate_rna_col='seurat_clusters')

cv_metrics <- cv_grn(test_srt, genes=genes_use, method='xgb', fit_intercept=T, alpha=1)
cv_metrics <- cv_grn(test_srt, genes=genes_use, method='bagging_ridge', fit_intercept=T, alpha=1)

test_srt <- infer_grn(test_srt, genes=genes_use, method='xgb', nrounds=100, nthread=-1)

test_srt <- infer_grn(test_srt, genes=genes_use, method='bagging_ridge', alpha=0.5, p_method='wilcox')
test_srt <- infer_grn(test_srt, genes=genes_use, method='bayesian_ridge', alpha=0.5, p_method='wilcox')

test_srt <- find_modules(test_srt, min_genes_per_module=0)

Params(test_srt)
NetworkParams(test_srt)

test_srt@grn@networks

coef(test_srt)
coef(test_srt, network='bayesian_ridge_network')
coef(test_srt, network='bagging_ridge_network')
gof(test_srt, network='bayesian_ridge_network')
gof(test_srt, network='bagging_ridge_network')
gof(test_srt)
modules <- NetworkModules(test_srt)

class(modules@meta)

aggregate_assay(test_srt, 'peaks_snn_res.50')


np <- 100
x <- seq(np) + rpois(np, 3)
z <- seq(np) + rpois(np, 5)
y <- rpois(np, 7) + x * 3 + z * 1.5

tbl <- tibble(
    x = x,
    y = y,
    z = z
)

strata <- sample(rep(LETTERS, np))[1:np]
flds <- cv_folds(tbl, strata=strata)

formula <- y ~ x + z
train <- tbl[-flds[[1]], ]
test <- tbl[flds[[1]], ]

score_glmnet(formula, train, test)
score_cvglmnet(formula, train, test)
score_bagging_ridge(formula, train, test)
cv_model(formula, tbl, method = 'bayesian_ridge', k_folds=5)

fit <- glmnetUtils::glmnet(formula, data=train)
y_true <- test[[formula[[2]]]]
y_pred <- predict(fit, newdata=test)[,1]

compute_metrics(y_true, y_pred)








