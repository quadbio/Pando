context('Tests for core Pando functions')

library(testthat)
library(readr)
library(GenomicRanges)
library(matrixStats)
library(doParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Pando)
library(Seurat)

data(motifs)

test_srt <- read_rds('../../../data/test_seurat.rds')
genes_use <- c('NEUROD6', 'POU5F1')

test_srt <- initiate_grn(test_srt, exclude_exons=F)
VariableFeatures(test_srt, assay='RNA') <- genes_use

test_that('initiate_grn returns SeuratPlus object.', {
    expect_equal(class(test_srt)[1], 'SeuratPlus')
})

test_that('initiate_grn creates RegualtoryNetwork object.', {
    expect_equal(class(GetGRN(test_srt))[1], 'RegulatoryNetwork')
})

test_srt <- find_motifs(
    test_srt, pfm=motifs[800:802], genome=BSgenome.Hsapiens.UCSC.hg38
)

test_that('find_motifs returns Motif object.', {
    expect_equal(class(NetworkRegions(test_srt)@motifs)[1], 'Motif')
})

test_srt <- infer_grn(test_srt, verbose=F)

test_that('infer_grn returns coefs.', {
    expect_setequal(unique(coef(test_srt)$target), c('NEUROD6', 'POU5F1'))
    expect_equal(dim(coef(test_srt))[1], 35)
})


test_that('infer_grn applies filtering.', {
    tst <- infer_grn(test_srt, verbose=F)
    expect_equal(dim(coef(tst))[1], 35)
    tst <- infer_grn(test_srt, verbose=F, tf_cor=0.2)
    expect_equal(dim(coef(tst))[1], 22)
    tst <- infer_grn(test_srt, verbose=F, tf_cor=0.9)
    expect_equal(dim(coef(tst))[1], 0)
    tst <- infer_grn(test_srt, verbose=F, peak_cor=0.2)
    expect_equal(dim(coef(tst))[1], 10)
    tst <- infer_grn(test_srt, verbose=F, peak_cor=0.9)
    expect_equal(dim(coef(tst))[1], 0)
})


test_that('initiate_grn excludes exons.', {
    tst <- initiate_grn(test_srt, exclude_exons=T)
    tst <- find_motifs(tst, pfm=motifs[800:802], genome=BSgenome.Hsapiens.UCSC.hg38)
    tst <- infer_grn(tst, verbose=F)
    expect_equal(dim(coef(tst))[1], 25)
})

test_that('getters return the expected type.', {
    expect_equal(class(GetNetwork(test_srt))[1], 'Network')
    expect_equal(class(NetworkFeatures(test_srt)), 'character')
    expect_equal(class(NetworkTFs(test_srt))[1], 'dgCMatrix')
    expect_equal(class(NetworkRegions(test_srt))[1], 'Regions')
    expect_equal(class(GetGRN(test_srt))[1], 'RegulatoryNetwork')
    expect_equal(class(NetworkModules(test_srt))[1], 'Modules')
    expect_equal(class(NetworkParams(test_srt)), 'list')
    expect_setequal(class(coef(test_srt)), c('tbl', 'tbl_df', 'data.frame'))
    expect_setequal(class(gof(test_srt)), c('tbl', 'tbl_df', 'data.frame'))
})


test_that('infer_grn makes new Network object and names correctly.', {
    tst <- infer_grn(test_srt, verbose=F, method='glmnet')
    expect_equal(length(GetGRN(tst)@networks), 2)
    expect_setequal(names(GetGRN(tst)@networks), c('glmnet_network', 'glm_network'))
})


registerDoParallel(3)

test_that('infer_grn with diverse parameters does not through errors.', {
    expect_error(infer_grn(test_srt, verbose=F, parallel=T), NA)
    expect_error(infer_grn(test_srt, verbose=F, parallel=T, tf_cor=0.2), NA)
    expect_error(infer_grn(test_srt, parallel=T, tf_cor=0.9), NA)
    expect_error(infer_grn(test_srt, verbose=F, parallel=T, peak_cor=0.1), NA)
    expect_error(infer_grn(test_srt, verbose=F, only_tss=T), NA)
    expect_error(infer_grn(test_srt, verbose=F, parallel=T, method='glmnet'), NA)
    expect_error(infer_grn(test_srt, verbose=F, method='cv.glmnet', alpha=1), NA)
    expect_error(infer_grn(test_srt, verbose=F, parallel=T, downstream=2000, upstream=0), NA)
    expect_error(infer_grn(test_srt, verbose=F, downstream=0, upstream=8000), NA)
    expect_error(infer_grn(test_srt, verbose=F, parallel=T, interaction_term='+'), NA)
    expect_error(infer_grn(test_srt, verbose=F, parallel=T, interaction_term='*'), NA)
    expect_error(infer_grn(test_srt, verbose=F, parallel=T, family='poisson'), NA)
    expect_error(infer_grn(test_srt, verbose=F, method='brms'), NA)
    expect_error(infer_grn(test_srt, verbose=F, method='xgb'), NA)
    expect_error(infer_grn(test_srt, verbose=F, method='bagging_ridge', p_method='wilcox'), NA)
    expect_error(infer_grn(test_srt, verbose=F, method='bagging_ridge', alpha=0.5), NA)
    expect_error(infer_grn(test_srt, verbose=F, method='bayesian_ridge'), NA)
})


test_that('infer_grn returns coefs for all models.', {
    tst <- infer_grn(test_srt, verbose=F, method='glm')
    expect_gt(dim(coef(tst, network='glm_network'))[1], 1)
    tst <- infer_grn(test_srt, verbose=F, method='glmnet')
    expect_gt(dim(coef(tst, network='glmnet_network'))[1], 1)
    tst <- infer_grn(test_srt, verbose=F, method='cv.glmnet')
    expect_gt(dim(coef(tst, network='cv.glmnet_network'))[1], 1)
    tst <- infer_grn(test_srt, verbose=F, method='brms')
    expect_gt(dim(coef(tst, network='brms_network'))[1], 1)
    tst <- infer_grn(test_srt, verbose=F, method='xgb')
    expect_gt(dim(coef(tst, network='xgb_network'))[1], 1)
    tst <- infer_grn(test_srt, verbose=F, method='bagging_ridge')
    expect_gt(dim(coef(tst, network='bagging_ridge_network'))[1], 1)
    tst <- infer_grn(test_srt, verbose=F, method='bayesian_ridge')
    expect_gt(dim(coef(tst, network='bayesian_ridge_network'))[1], 1)
})


test_that('find_modules works for glm models.', {
    tst <- infer_grn(test_srt, verbose=F, parallel=T)
    expect_error(find_modules(tst), NA)
    tst <- find_modules(tst)
    expect_true('tbl'%in%class(NetworkModules(tst)@meta))
})

test_that('find_modules works for glmnet models.', {
    tst <- infer_grn(test_srt, verbose=F, parallel=T, method='cv.glmnet')
    expect_error(find_modules(tst), NA)
    tst <- find_modules(tst)
    expect_true('tbl'%in%class(NetworkModules(tst)@meta))
    tst <- infer_grn(test_srt, verbose=F, parallel=T, method='glmnet')
    expect_error(find_modules(tst), NA)
    tst <- find_modules(tst)
    expect_true('tbl'%in%class(NetworkModules(tst)@meta))
})

test_that('find_modules works for brms models.', {
    tst <- infer_grn(test_srt, verbose=F, parallel=T, method='brms')
    expect_error(find_modules(tst), NA)
    tst <- find_modules(tst)
    expect_true('tbl'%in%class(NetworkModules(tst)@meta))
})


