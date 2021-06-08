context('Tests for core Pando functions')

library(testthat)
library(readr)
library(GenomicRanges)
library(matrixStats)
library(doParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Pando)

data(motifs)

test_srt <- read_rds('../../../data/test_seurat.rds')
test_srt <- initiate_grn(test_srt, genes = c('NEUROD6', 'POU5F1', 'PAX6'))

test_that('initiate_grn returns SeuratPlus object.', {
    expect_equal(class(test_srt)[1], 'SeuratPlus')
})

test_that('initiate_grn creates RegualtoryNetwork object.', {
    expect_equal(class(GetNetworkData(test_srt))[1], 'RegulatoryNetwork')
})

test_that('initiate_grn filters target genes.', {
    expect_setequal(NetworkFeatures(test_srt), c('NEUROD6', 'POU5F1'))
})


test_srt <- find_motifs(
    test_srt, pfm=motifs[800:802], genome=BSgenome.Hsapiens.UCSC.hg38
)

test_that('find_motifs returns Motif object.', {
    expect_equal(class(NetworkRegions(test_srt)@motifs)[1], 'Motif')
})


test_srt <- infer_grn(test_srt, verbose=F)

test_that('infer_grn returns Network object.', {
    expect_equal(class(GetNetworkData(test_srt)@network)[1], 'Network')
})


test_that('infer_grn returns coefs.', {
    expect_setequal(unique(coef(test_srt)$target), c('NEUROD6', 'POU5F1'))
    expect_equal(dim(coef(test_srt)), c(37, 6))
})


test_that('infer_grn applies filtering.', {
    tst <- infer_grn(test_srt, verbose=F)
    expect_equal(dim(coef(tst)), c(37, 6))
    tst <- infer_grn(test_srt, verbose=F, tf_cor=0.2)
    expect_equal(dim(coef(tst)), c(24, 6))
    tst <- infer_grn(test_srt, verbose=F, tf_cor=0.9)
    expect_equal(dim(coef(tst)), c(0, 0))
    tst <- infer_grn(test_srt, verbose=F, peak_cor=0.2)
    expect_equal(dim(coef(tst)), c(11, 6))
    tst <- infer_grn(test_srt, verbose=F, peak_cor=0.9)
    expect_equal(dim(coef(tst)), c(0, 0))
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
})





