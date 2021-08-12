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

registerDoParallel(1)

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
test_srt <- infer_grn(test_srt, peak_to_gene_method = 'GREAT', parallel=T, method = 'cv.glmnet', nlambda=100)
test_srt <- infer_grn(test_srt, peak_to_gene_method = 'GREAT', method = 'brms', backend='cmdstanr', prior=horseshoe())

test_srt <- find_modules(test_srt, min_genes_per_module=0)

coef(test_srt)
NetworkModules(test_srt)



#### Try bayesian regression ####
np <- 100
x <- seq(np) + rnorm(np, sd=3)
z <- seq(np) + rnorm(np, sd=5)
y <- round(rnorm(np, sd=50, mean=1000) + x * 3 + z * 1.5)

tbl <- tibble(
    x = x,
    y = y,
    z = z
)


p <- prior(normal(0,5))
# p <- prior(horseshoe(scale_global=1000))

brm_fit <- suppressMessages(brm(y ~ x : z, data=tbl, backend='cmdstanr', family=gaussian, prior=p, silent=T, refresh = 0))

class(brm_fit)


fixef(brm_fit)
brm_smry <- summary(brm_fit, )
bayes_R2(brm_fit)
pvals <- p_map(brm_fit)$p_MAP







