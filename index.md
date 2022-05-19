# Pando <img src="man/figures/logo.png" align="right" width="180"/>

Pando leverages multi-modal singel-cell measurements to infer gene regulatory networks (GRNs) using a flexible modelling framework. By modeling the relationship between TF-binding site pairs with the expression of target genes, Pando simultaneously infers gene modules and sets of regulatory regions for each transcription factor.


## Installation

```r
devtools::install_github('quadbiolab/Pando')
```


## About Pando

The fate and state of a cell are regulated through complex circuits of transcription factors (TFs) converging at regulatory elements to enable precise control of gene expression. Modern single-cell genomic approaches allow the simultaneous profiling of gene expression and chromatin accessibility in individual cells, which opens up new opportunities for the inference of gene regulatory networks (GRNs). 

The unifying idea behind many modern GRN inference methods is to model the expression of each gene as a function of TF abundances. The weights or coefficients of this model can then be interpreted as a measure of the regulatory interaction between TF and target gene. Additional (epi-) genomic information (such as predicted TF binding sites) is often used to constrain or refine the model.

Pando tries to generalize this concept to make use of the multi-modal nature of modern single-cell technologies by incorporating TF binding information directly into the model. By utilizing jointly measured or integrated scRNA-seq and scATAC-seq data, Pando models the expression of genes based on the interaction of TF expression with the accessibility of their putative binding site. By offering a number of different pre-processing and modeling choices, Pando strives to be a modular and flexible framework for single-cell GRN inference.



## Usage overview

### Initiating the GRN
Pando interacts directly with [Seurat objects](https://satijalab.org/seurat/) and integrates well with [Seurat](https://satijalab.org/seurat/) and [Signac](https://satijalab.org/signac/) workflows. To use Pando, you'll need a Seurat object with two assays, one with scRNA-seq transcript counts and one with scATAC-seq peak accessibility. With this object (let's call it `seurat_object`) ready, you can start off by initializing the GRN using the function `initiate_grn()`:

```r
seurat_object <- initiate_grn(seurat_object)
```

This will create a `RegulatoryNetwork` object inside the Seurat object and select candidate regulatory regions. By default, Pando will consider all peaks as putative regulatory regions, but the set of candidate regions can be constrained by providing a `GenomicRanges` object in the `regions` argument. Pando ships with a set of conserved regions (`phastConsElements20Mammals.UCSC.hg38`) as well as predicted regulatory elements from [ENCODE](https://screen.encodeproject.org/) (`SCREEN.ccRE.UCSC.hg38`) for the human genome (hg38), which could be used here. However, one could also select candidate regions in other ways, for instance by using [Cicero](https://cole-trapnell-lab.github.io/cicero-release/).

### Finding TF binding sites
Once the `RegulatoryNetwork` object is initiated with candidate regions, we can scan for TF binding motifs in these regions by using the function `find_motifs()`

```r
library(BSgenome.Hsapiens.UCSC.hg38)
data(motifs)

seurat_object <- find_motifs(
    seurat_object,
    pfm = motifs,
    genome = BSgenome.Hsapiens.UCSC.hg38
)
```

This uses [motifmatchr](https://github.com/GreenleafLab/motifmatchr) to pair up TFs with their putative binding sites. Pando provides a custom motif database (`motifs`) compiled from [JASPAR](http://jaspar.genereg.net/) and [CIS-BP](http://cisbp.ccbr.utoronto.ca/), but in principle any `PFMatrixList` object can be provided here. A data frame with motif-to-TF assignments can be provided in the `motif_tfs` argument.

### Inferring the GRN
Now everything should be ready to infer the GRN by fitting regression models for the expression of each gene. In Pando, this can be done by using the function `infer_grn()`:

```r
seurat_object <- infer_grn(
    seurat_object,
    peak_to_gene_method = 'Signac',
    method = 'glm'
)
```

Here, we first select regions near genes, either by simply considering a distance upstream and/or downstream of the gene (`peak_to_gene_method='Signac'`) or by also considering overlapping regulatory regions as is done by [GREAT](http://great.stanford.edu/public/html/) (`peak_to_gene_method='GREAT'`). 

You can also choose between a number of different models using the `method` argument, such as GLMs (`'glm'`) regularized GLMs [(`'glmnet'`, `'cv.glmnet'`)](https://glmnet.stanford.edu/articles/glmnet.html) or Bayesian regression models [(`'brms'`)](https://paul-buerkner.github.io/brms/). We are also integrated gradient boosting regression with [XGBoost](https://xgboost.readthedocs.io/en/latest/R-package/xgboostPresentation.html) as it is used by [GRNBoost](https://github.com/aertslab/GRNBoost)/[SCENIC](https://scenic.aertslab.org/) as well as bagging and bayesian ridge models with [scikit-learn](https://scikit-learn.org/stable/) as they are used by [CellOracle](https://github.com/morris-lab/CellOracle).

Once the models are fit, model coefficients can be inspected with

```r
coef(seurat_object)
```


### Module extraction
Based on the model coefficients, we can construct a network between TFs and target genes. This can be further summarized to construct gene and regulatory modules with the set of target genes and regulatory regions for each TF. In Pando we do this with 

```r
seurat_object <- find_modules(seurat_object)
```

To access the extracted modules, you can use the function `NetworkModules()`:

```
modules <- NetworkModules(seurat_object)
modules@meta
```

The `meta` slot holds a dataframe with module inforamtion. 

If you are curious to find out more, check out our vignettes!










