# Pando <img src="man/figures/logo.png" align="right" width="180"/>

Pando leverages multi-modal singel-cell measurements to infer gene regulatory networks (GRNs) using a flexible modelling framework. By modeling the relationship between TF-binding site pairs with the expression of target genes, Pando simultaneously infers gene modules and sets of regulatory regions for each transcription factor.


## Installation

```r
devtools::install_github('quadbiolab/Pando')
```


## About Pando

The fate and state of a cell is regulated through complex circuits of transcription factors (TFs) converging at regulatory elements to enable precise control of gene expression. Modern single-cell genomic approaches allow the simultaneous profiling of gene expression and chromatin accessibility in individual cells, which opens up new opportunities for the inference of gene regulatory networks (GRNs). 

The unifying idea behind many modern GRN inference methods is to model the expression of each gene as a function of TF abundances. The weights or coefficients of this model can then be interpreted as a measure of the regulatory interaction between TF and target gene. Additional (epi-) genomic information (such as predicted TF binding sites) is often used to constrain or refine the model.

Pando tries to generalize this concept to make use of the multi-modal nature of modern single-cell technologies by incorporating TF binding information directly into the model. By utilizing jointly measured or integrated scRNA-seq and scATAC-seq data, Pando models the expression of genes based the interaction of TF expression with the accessibility of their putative binding site. By offering a number different pre-processing and modelling choices, Pando strives to be a modular and flexible framework for single-cell GRN inference.



## Usage overview

Pando interacts directly with [Seurat objects](https://satijalab.org/seurat/) and integrates well [Seurat](https://satijalab.org/seurat/) and [Signac](https://satijalab.org/signac/) workflows. To use Pando, you'll need a Seurat object with two assays, one scRNA-seq transcript counts and one with scATAC-seq peak accessibility. With this object (let's call it `seurat_object`) ready, you can start off by inizializing the GRN using the function `initiate_grn()`:

```r
seurat_object <- initiate_grn(seurat_object)
```

This will create a `RegulatoryNetwork` object inside the Seurat object and select candidate regulatory regions. Per default, Pando will consider all peaks as putative regulatory regions, but the set of candidate regions can be constrained by providing a `GenomicRanges` object in the `regions` argument. 


