library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Pando)
library(doParallel)

data(motifs)
data(motif2tf)

registerDoParallel(4)

test_srt <- read_rds('../data/test_seurat.rds')
test_srt <- Seurat::FindVariableFeatures(test_srt, assay='RNA')
test_srt <- initiate_grn(test_srt)
test_srt <- find_motifs(
    test_srt,
    pfm = motifs[1:10],
    genome = BSgenome.Hsapiens.UCSC.hg38
)
test_srt <- infer_grn(test_srt, parallel = T)
coef(test_srt)


Signac::Annotation(test_srt)

test_srt@network@regions@motifs
names(motifs)

ll <- list(
    a = 1,
    b= 3,
    c=4,
    d=9
)

vv <- c(1,2,3,4,9,8)
names(vv) <- LETTERS[1:6]

flexapply(ll, sqrt, parallel=F)

mm <- matrix(1:100, nrow=10)
rownames(mm) <- LETTERS[1:10]
colnames(mm) <- letters[1:10]

mmr <- presto::rank_matrix(mm)$X_ranked
sparse_cor(mmr)


