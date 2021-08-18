#' @import dplyr tibble
NULL


#' Infer a Gene Regulatory Network with \code{Pando}
#'
#' @import Matrix
#' @import sparseMatrixStats
#' @importFrom purrr map_dfr map_lgl map_dbl map
#' @importFrom stringr str_replace_all
#'
#' @param genes A character vector with the target genes to consider for GRN inference.
#' Takes all VariableFeatures in the object per default.
#' @param peak_to_gene_method Character specifying the method to
#' link peak overlapping motif regions to nearby genes. One of 'Signac or 'GREAT'.
#' @param upstream Integer defining the distance upstream of the gene to consider as potential regulatory region.
#' @param downstream Integer defining the distance downstream of the gene to consider as potential regulatory region.
#' @param extend Integer defining the distance from the upstream and downstream of the basal regulatory region.
#' Only used of `peak_to_gene_method = 'GREAT'`.
#' @param only_tss Logical. Measure distance from the TSS (\code{TRUE}) or from the entire gene body (\code{FALSE}).
#' @param parallel Logical. Whether to parellelize the computation with \code{\link[foreach]{foreach}}.
#' @param tf_cor Threshold for TF - target gene correlation.
#' @param peak_cor Threshold for binding peak - target gene correlation.
#' @param method A character string indicating the method to fit the model.
#' Possible values are \code{'glm'}, \code{'glmnet'} and \code{'cv.glmnet'}.
#' @param alpha The elasticnet mixing parameter. See \code{\link[glmnet]{glmnet}} for details.
#' @param family A description of the error distribution and link function to be used in the model.
#' See \code{\link[family]{stats}} for mode details.
#' @param interaction_term The interaction term to use in the model between TF and binding site.
#' * \code{'+'} for additive interaction.
#' * \code{':'} for 'multiplicative' interaction.
#' * \code{'\*'} for crossing interaction, i.e. additive AND 'multiplicative'.
#' For more info, see \code{\link[formula]{stats}}
#' @param adjust_method Method fro adjusting p-values.
#' @param verbose Logical. Display messages
#' @param ... Other parameters for the model fitting function.
#'
#' @return A SeuratPlus object.
#'
#' @rdname infer_grn
#' @export
#' @method infer_grn SeuratPlus
infer_grn.SeuratPlus <- function(
    object,
    genes = NULL,
    network_name = paste0(method, '_network'),
    peak_to_gene_method = 'Signac',
    upstream = 100000,
    downstream = 0,
    extend = 1000000,
    only_tss = FALSE,
    parallel = FALSE,
    tf_cor = 0.1,
    peak_cor = 0.,
    method = 'glm',
    alpha = 0.5,
    family = 'gaussian',
    interaction_term = ':',
    adjust_method = 'fdr',
    verbose = TRUE,
    ...
){
    # Get variables from object
    params <- Params(object)
    motif2tf <- NetworkTFs(object)
    if (is.null(motif2tf)){
        stop('Motif matches have not been found. Please run find_motifs() first.')
    }
    gene_annot <- Signac::Annotation(object[[params$peak_assay]])
    if (is.null(gene_annot)){
        stop('Please provide a gene annotation for the ChromatinAssay.')
    }
    # Select target genes for GRN inference
    if (is.null(genes)){
        genes <- VariableFeatures(object, assay=params$rna_assay)
        if (is.null(genes)){
            stop('Please provide a set of features or run FindVariableFeatures()')
        }
    }
    # Select genes to use by intersecting annotated genes with all
    # detected genes in the object
    features <- intersect(gene_annot$gene_name, genes) %>%
        intersect(rownames(GetAssay(object, params$rna_assay)))
    gene_annot <- gene_annot[gene_annot$gene_name%in%features, ]

    # Get regions
    regions <- NetworkRegions(object)
    log_message('Selecting candidate regulatory regions near genes', verbose=verbose)
    peak_data <- t(Seurat::GetAssayData(object, assay=params$peak_assay)[regions@peaks, ])
    colnames(peak_data) <- rownames(regions@motifs@data)
    peaks2motif <- regions@motifs@data

    # Find candidate regions near gene bodies
    peaks_near_gene <- find_peaks_near_genes(
        peaks = regions@ranges,
        method = peak_to_gene_method,
        genes = gene_annot,
        upstream = upstream,
        downstream = downstream,
        only_tss = only_tss
    )
    peaks2gene <- aggregate_matrix(t(peaks_near_gene), groups=colnames(peaks_near_gene), fun='sum')

    # Select peaks passing criteria
    peaks_at_gene <- as.logical(colMaxs(peaks2gene))
    peaks_with_motif <- as.logical(rowMaxs(peaks2motif))

    # Subset data to good peaks
    peaks_use <- peaks_at_gene & peaks_with_motif
    peaks2gene <- peaks2gene[, peaks_use]
    peaks2motif <- peaks2motif[peaks_use, ]
    peak_data <- peak_data[, peaks_use]

    log_message('Preparing model input', verbose=verbose)
    gene_data <- t(Seurat::GetAssayData(object, assay=params$rna_assay))
    tfs_use <- colnames(motif2tf)
    motif2tf <- motif2tf[, tfs_use]

    log_message('Fitting models for ', length(features), ' target genes' , verbose=verbose)
    names(features) <- features
    model_fits <- map_par(features, function(g){

        # Select peaks near gene
        gene_peaks <- as.logical(peaks2gene[g, ])
        if (sum(gene_peaks)==0){
            log_message('Warning: No peaks found near ', g, verbose=verbose==2)
            return()
        }

        # Select peaks correlating with target gene expression
        peak_x <- peak_data[, gene_peaks, drop=F]
        peak_g_cor <- sparse_cor(peak_x, gene_data[, g, drop=F])
        peak_g_cor[is.na(peak_g_cor)] <- 0
        peaks_use <- rownames(peak_g_cor)[abs(peak_g_cor[, 1]) > peak_cor]
        if (length(peaks_use)==0){
            log_message('Warning: No correlating peaks found for ', g, verbose=verbose==2)
            return()
        }
        peak_x <- peak_x[, peaks_use, drop=F]
        peak_motifs <- peaks2motif[gene_peaks, , drop=F][peaks_use, , drop=F]

        # Select TFs with motifs in peaks
        gene_peak_tfs <- map(rownames(peak_motifs), function(p){
            x <- as.logical(peak_motifs[p, ])
            peak_tfs <- colMaxs(motif2tf[x, , drop=F])
            peak_tfs <- colnames(motif2tf)[as.logical(peak_tfs)]
            peak_tfs <- setdiff(peak_tfs, g)
            return(peak_tfs)
        })
        names(gene_peak_tfs) <- rownames(peak_motifs)

        # Check correlation of peaks with target gene
        gene_tfs <- purrr::reduce(gene_peak_tfs, union)
        tf_x <- gene_data[, gene_tfs, drop=F]
        tf_g_cor <- sparse_cor(tf_x, gene_data[, g, drop=F])
        tf_g_cor[is.na(tf_g_cor)] <- 0
        tfs_use <- rownames(tf_g_cor)[abs(tf_g_cor[, 1]) > tf_cor]
        if (length(tfs_use)==0){
            log_message('Warning: No correlating TFs found for ', g, verbose=verbose==2)
            return()
        }

        # Filter TFs and make formula string
        frml_string <- map(names(gene_peak_tfs), function(p){
            peak_tfs <- gene_peak_tfs[[p]]
            peak_tfs <- peak_tfs[peak_tfs%in%tfs_use]
            if (length(peak_tfs)==0){
                return()
            }
            peak_name <- str_replace_all(p, '-', '_')
            tf_name <- str_replace_all(peak_tfs, '-', '_')
            formula_str <- paste(
                paste(peak_name, interaction_term, tf_name, sep=' '), collapse = ' + ')
            return(list(tfs=peak_tfs, frml=formula_str))
        })
        frml_string <- frml_string[!map_lgl(frml_string, is.null)]
        if (length(frml_string)==0){
            log_message('Warning: No valid peak:TF pairs found for ', g, verbose=verbose==2)
            return()
        }

        target <- str_replace_all(g, '-', '_')
        model_frml <- as.formula(
            paste0(target, ' ~ ', paste0(map(frml_string, function(x) x$frml),  collapse=' + '))
        )

        # Get expression data
        nfeats <- sum(map_dbl(frml_string, function(x) length(x$tfs)))
        gene_tfs <- purrr::reduce(map(frml_string, function(x) x$tfs), union)
        gene_x <- gene_data[, union(g, gene_tfs), drop=F]
        model_mat <- as.data.frame(cbind(gene_x, peak_x))
        colnames(model_mat) <- str_replace_all(colnames(model_mat), '-', '_')

        log_message('Fitting model with ', nfeats, ' variables for ', g, verbose=verbose==2)
        fit <- try(fit_model(
            model_frml,
            data = model_mat,
            family = family,
            method = method,
            alpha = alpha,
            ...
        ), silent=TRUE)
        if (any(class(fit)=='try-error')){
            log_message('Warning: Fitting model failed for ', g, verbose=verbose)
            log_message(fit, verbose=verbose==2)
            return()
        } else {
            fit$gof$nvariables <- nfeats
            return(fit)
        }
    }, verbose=verbose, parallel=parallel)
    model_fits <- model_fits[!map_lgl(model_fits, is.null)]
    if (length(model_fits)==0){
        log_message('Warning: Fitting model failed for all genes.', verbose=verbose)
    }
    coefs <- map_dfr(model_fits, function(x) x$coefs, .id='target')
    coefs <- format_coefs(coefs, term=interaction_term, adjust_method=adjust_method)
    gof <- map_dfr(model_fits, function(x) x$gof, .id='target')

    params <- list()
    params[['method']] <- method
    params[['family']] <- family
    params[['dist']] <- c('upstream'=upstream, 'downstream'=downstream)
    params[['only_tss']] <- only_tss
    params[['interaction']] <- interaction_term
    params[['tf_cor']] <- tf_cor
    params[['peak_cor']] <- peak_cor

    network_obj <- new(
        Class = 'Network',
        features = features,
        coefs = coefs,
        fit = gof,
        params = params
    )
    object@grn@networks[[network_name]] <- network_obj
    return(object)
}



#' Format network coefficients
#'
#' @import stringr
#'
#' @param coefs A data frame with coefficients
#'
#' @return A data frame.
#'
#' @export
format_coefs <- function(coefs, term=':', adjust_method='fdr'){

    if (dim(coefs)[1] == 0){
        return(coefs)
    }

    if ('pval' %in% colnames(coefs)){
        coefs$padj <- p.adjust(coefs$pval, method=adjust_method)
    }

    term_pattern <- paste0('(.+)', term, '(.+)')
    region_pattern <- '[\\d\\w]+_\\d+_\\d+'
    coefs_use <- coefs %>%
        filter(!term%in%c('(Intercept)', 'Intercept')) %>%
        mutate(
            tf_ = str_replace(term, term_pattern, '\\1'),
            region_ = str_replace(term, term_pattern, '\\2')
        ) %>%
        mutate(
            tf = ifelse(str_detect(tf_, region_pattern), region_, tf_),
            region = ifelse(!str_detect(tf_, region_pattern), region_, tf_)
        ) %>%
        select(-region_, -tf_) %>%
        mutate(
            region = str_replace_all(region, '_', '-'),
            tf = str_replace_all(tf, '_', '-'),
            target = str_replace_all(target, '_', '-')
        ) %>%
        select(tf, target, region, term, everything())
    return(coefs_use)
}



#' Find TF modules in regulatory network
#'
#' @importFrom purrr map map_chr
#' @importFrom stringr str_split str_replace_all
#'
#' @param p_thresh Float indicating the significance threshold on the adjusted p-value.
#' @param rsq_thresh Float indicating the \eqn{R^2} threshold on the adjusted p-value.
#' @param nvar_thresh Integer indicating the minimum number of variables in the model.
#' @param min_genes_per_module Integer indicating the minimum number of genes in a module.
#'
#' @return A RegulatoryNetwork object.
#'
#' @rdname find_modules
#' @export
#' @method find_modules RegulatoryNetwork
find_modules.RegulatoryNetwork <- function(
    object,
    p_thresh = 0.05,
    rsq_thresh = 0.1,
    nvar_thresh = 10,
    min_genes_per_module = 5
){
    fit_method <- Params(object)$method

    models_use <- gof(object) %>%
        filter(rsq>rsq_thresh & nvariables>nvar_thresh) %>%
        pull(target) %>%
        unique()

    modules <- coef(object) %>%
        filter(target %in% models_use)

    if (fit_method %in% c('cv.glmnet', 'glmnet')){
        modules <- modules %>%
            filter(estimate != 0)
    } else {
        modules <- modules %>%
            filter(ifelse(is.na(padj), T, padj<p_thresh))
    }

    modules <- modules %>%
        group_by(target) %>%
        mutate(nvars=n()) %>%
        group_by(target, tf) %>%
        mutate(tf_sites_per_gene=n()) %>%
        group_by(target) %>%
        mutate(
            tf_per_gene=length(unique(tf)),
            peak_per_gene=length(unique(region))
        ) %>%
        group_by(tf) %>%
        mutate(gene_per_tf=length(unique(target))) %>%
        group_by(target, tf)

    if (fit_method %in% c('cv.glmnet', 'glmnet')){
        modules <- modules %>%
            summarize(
                estimate=sum(estimate),
                n_regions=peak_per_gene,
                n_genes=gene_per_tf,
                n_tfs=tf_per_gene,
                regions=paste(region, collapse=';')
            )
    } else {
        modules <- modules %>%
            summarize(
                estimate=sum(estimate),
                n_regions=peak_per_gene,
                n_genes=gene_per_tf,
                n_tfs=tf_per_gene,
                regions=paste(region, collapse=';'),
                pval=min(pval),
                padj=min(padj)
            )
    }

    modules <- modules %>%
        distinct() %>%
        arrange(tf)

    module_pos <- modules %>%
        filter(estimate>0) %>%
        group_by(tf) %>% filter(n()>min_genes_per_module) %>%
        group_split() %>% {names(.) <- map_chr(., function(x) x$tf[[1]]); .} %>%
        map(function(x) x$target)

    module_neg <- modules %>%
        filter(estimate<0) %>%
        group_by(tf) %>% filter(n()>min_genes_per_module) %>%
        group_split() %>% {names(.) <- map_chr(., function(x) x$tf[[1]]); .} %>%
        map(function(x) x$target)

    regions_pos <- modules %>%
        filter(estimate>0) %>%
        group_by(tf) %>% filter(n()>min_genes_per_module) %>%
        group_split() %>% {names(.) <- map_chr(., function(x) x$tf[[1]]); .} %>%
        map(function(x) unlist(str_split(x$regions, ';')))

    regions_neg <- modules %>%
        filter(estimate<0) %>%
        group_by(tf) %>% filter(n()>min_genes_per_module) %>%
        group_split() %>% {names(.) <- map_chr(., function(x) x$tf[[1]]); .} %>%
        map(function(x) unlist(str_split(x$regions, ';')))

    module_feats <- list(
        'genes_pos' = module_pos,
        'genes_neg' = module_neg,
        'regions_pos' = regions_pos,
        'regions_neg' = regions_neg
    )

    object@network@modules@meta <- modules
    object@network@modules@features <- module_feats

    return(object)
}



#' @importFrom purrr map
#'
#' @return A Network object.
#'
#' @rdname find_modules
#' @export
#' @method find_modules SeuratPlus
find_modules.SeuratPlus <- function(
    object,
    network = 'glm_network',
    p_thresh = 0.05,
    rsq_thresh = 0.1,
    nvar_thresh = 10,
    min_genes_per_module = 5
){
    params <- Params(object)
    regions <- NetworkRegions(object)
    network <- GetNetwork(object, network=network)
    network <- find_modules(
        network,
        p_thresh = p_thresh,
        rsq_thresh = rsq_thresh,
        nvar_thresh = nvar_thresh,
        min_genes_per_module = min_genes_per_module
    )
    modules <- NetworkModules(network)

    reg2peaks <- rownames(GetAssay(object, assay=params$peak_assay))[regions@peaks]
    names(reg2peaks) <- Signac::GRangesToString(regions@ranges)
    peaks_pos <- modules@features$regions_pos %>% map(function(x) unique(reg2peaks[x]))
    peaks_neg <- modules@features$regions_neg %>% map(function(x) unique(reg2peaks[x]))
    modules@features[['peaks_pos']] <- peaks_pos
    modules@features[['peaks_neg']] <- peaks_neg

    object@grn@network@modules <- modules
    return(object)
}


