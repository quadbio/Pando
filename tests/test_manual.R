library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Pando)
library(doParallel)
library(devtools)
library(Signac)
library(Seurat)

data(motifs)
data(phastConsElements20Mammals.UCSC.hg38)
data(EnsDb.Hsapiens.v93.annot.UCSC.hg38)
data(motif2tf)

genes_use <- c('PAX6', 'POU5F1', 'DLX6', 'NFIA')
tfs_use <- c('PAX6', 'POU5F1', 'DLX6', 'NFIA')

motif2tf_use <- motif2tf %>%
    filter(origin=='JASPAR2020') %>%
    rename('mot'=motif)

registerDoParallel(4)

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
test_srt <- infer_grn(test_srt, peak_to_gene_method = 'GREAT', parallel=T)
test_srt <- find_modules(test_srt, min_genes_per_module=0)

coef(test_srt)
NetworkModules(test_srt)


















