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

registerDoParallel(2)

genes_use <- c('PAX6', 'POU5F1', 'DLX6', 'NFIA', 'NEUROD6', 'NFIB', 'RAX', 'FGF8', 'CRABP1', 'JUNB',
    'HOPX', 'LHX9', 'NEUROG2', 'LHX1', 'SIX6', 'FOXC1', 'DKK1', 'PCP4', 'SPRY1')
tfs_use <- unique(c(genes_use, 'PAX6', 'POU5F1', 'DLX6', 'NFIA', 'GLI3', 'NFIX', 'SOX3', 'SIX3'))

motif2tf_use <- motif2tf %>%
    filter(tf%in%tfs_use) %>%
    rename('mot'=motif)

test_srt <- read_rds('../data/test_seurat.rds')

test_srt <- initiate_grn(
    test_srt,
    regions = phastConsElements20Mammals.UCSC.hg38,
    exclude_exons = F
)

test_srt <- find_motifs(
    test_srt,
    motif_tfs = motif2tf_use,
    pfm = motifs[unique(motif2tf_use$mot)],
    genome = BSgenome.Hsapiens.UCSC.hg38
)


test_srt <- infer_grn(test_srt, genes=genes_use,
    peak_to_gene_method = 'GREAT', parallel=F)

test_srt <- find_modules(test_srt, min_genes_per_module=0, nvar_thresh=2)

test_srt <- get_network_graph(test_srt, n_neighbors=2)
plot_network_graph(test_srt, layout='umap')

NetworkGraph(test_srt)

plot_gof(test_srt, point_size=2)
plot_module_metrics(test_srt)

Params(test_srt)
NetworkParams(test_srt)
NetworkModules(test_srt)@meta
NetworkModules(test_srt)@params

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








