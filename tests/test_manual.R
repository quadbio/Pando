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
data(SCREEN.ccRE.UCSC.hg38)
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
test_srt <- infer_grn(test_srt, peak_to_gene_method = 'GREAT', method='xgb')

test_srt <- find_modules(test_srt, min_genes_per_module=0)

NetworkParams(test_srt)$method

coef(test_srt)
gof(test_srt)
modules <- NetworkModules(test_srt)

class(modules@meta)
















