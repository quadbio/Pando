context('Tests for Pando plotting functions')

library(BSgenome.Hsapiens.UCSC.hg38)
library(Pando)
library(Seurat)

data(motifs)

test_srt <- read_rds('../../../data/test_seurat.rds')
genes_use <- c('NEUROD6', 'POU5F1')

test_srt <- initiate_grn(test_srt, exclude_exons=F)
VariableFeatures(test_srt, assay='RNA') <- genes_use

test_srt <- find_motifs(
    test_srt, pfm=motifs[800:805], genome=BSgenome.Hsapiens.UCSC.hg38
)

test_srt <- infer_grn(test_srt, verbose=F)



