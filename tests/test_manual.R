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
    motif_tfs = motif2tf,
    pfm = motifs[unique(motif2tf_use$mot)],
    genome = BSgenome.Hsapiens.UCSC.hg38
)

test_srt <- infer_grn(test_srt, genes=genes_use, method='xgb',
    peak_to_gene_method = 'GREAT', parallel=F, verbose = 2)

annot <- Annotation(test_srt)

plot_module_metrics(test_srt)
plot_gof(test_srt)

annot <- annot[annot$gene_name%in%c('NEUROD6', 'MSX2', 'WLS', 'HES1')] %>%
    CollapseToLongestTranscript() %>%
    Extend(upstream = 200000, downstream = 200000)
annot$region_type <- 'TAD'
annot %>% write_rds('../data/example_tads.rds')

test_srt <- infer_grn(test_srt, genes=genes_use, peak_to_gene_domains=Annotation(test_srt))

test_srt <- find_modules(test_srt, min_genes_per_module=0, nvar_thresh=2)

test_srt <- get_network_graph(test_srt, n_neighbors=2, umap_method = 'weighted')
plot_network_graph(test_srt, layout='fr')

test_srt <- get_tf_network(test_srt, tf='NFIB')
plot_tf_network(test_srt, tf='NFIB', circular=F, label_nodes = 'tfs', edge_width = 3)

NetworkGraph(test_srt)

plot_gof(test_srt, point_size=2)
plot_module_metrics(test_srt)








