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

test_srt <- initiate_grn(
    test_srt,
    regions = phastConsElements20Mammals.UCSC.hg38
)
test_srt <- find_motifs(
    test_srt,
    motif_tfs = motif2tf_use,
    pfm = motifs[1:10],
    genome = BSgenome.Hsapiens.UCSC.hg38
)

test_srt <- aggregate_assay(test_srt, group_name = 'seurat_clusters')

test_srt <- infer_grn(test_srt, genes=genes_use,
    peak_to_gene_method = 'GREAT', parallel=F,
    aggregate_peaks_col='seurat_clusters', aggregate_rna_col='seurat_clusters')

test_srt <- infer_grn(test_srt, genes=genes_use, method='bayesian_ridge')
test_srt <- infer_grn(test_srt, genes=genes_use, method='bagging_ridge', n_jobs=1, p_method='wilcox')

test_srt <- find_modules(test_srt, min_genes_per_module=0)

Params(test_srt)
NetworkParams(test_srt)

test_srt@grn@networks

coef(test_srt)
coef(test_srt, network='bayesian_ridge_network')
coef(test_srt, network='bagging_ridge_network')
gof(test_srt)
modules <- NetworkModules(test_srt)

class(modules@meta)

aggregate_assay(test_srt, 'peaks_snn_res.50')


np <- 10000
x <- seq(np) + rpois(np, 3)
z <- seq(np) + rpois(np, 5)
y <- rpois(np, 7) + x * 3 + z * 1.5

tbl <- tibble(
    x = x,
    y = y,
    z = z
)

fit_bayesian_ridge(y ~ x + z, tbl)










