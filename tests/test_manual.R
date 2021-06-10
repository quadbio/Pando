library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Pando)
library(doParallel)
library(devtools)

data(motifs)
data(phastConsElements20Mammals.UCSC.hg38)
data(motif2tf)

registerDoParallel(4)

test_srt <- read_rds('../data/test_seurat.rds')
test_srt <- initiate_grn(
    test_srt,
    genes = c('NEUROD6', 'POU5F1', 'DLX6', 'NFIA'),
    regions = phastConsElements20Mammals.UCSC.hg38
)
test_srt <- find_motifs(
    test_srt,
    pfm = motifs[1:10],
    genome = BSgenome.Hsapiens.UCSC.hg38
)
test_srt <- infer_grn(test_srt, parallel=T)
test_srt <- find_modules(test_srt, min_genes_per_module = 0)

NetworkModules(test_srt)@meta
test_srt@grn@regions
format_coefs(coef(test_srt))


