#' @import dplyr tibble
NULL


#' Infer a Gene Regulatory Network with \code{Pando}
#'
#' @import Matrix
#' @import sparseMatrixStats
#' @importFrom purrr map_dfr map_lgl map_dbl map
#' @importFrom stringr str_replace_all
#'
#' @param upstream Integer defining the distance upstream of the gene to consider as potential regulatory region.
#' @param downstream Integer defining the distance downstream of the gene to consider as potential regulatory region.
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
#' @param verbose Logical. Display messages
#'
#' @return A SeuratPlus object.
#'
#' @rdname infer_grn
#' @export
#' @method infer_grn SeuratPlus
infer_grn.SeuratPlus <- function(
    object,
    upstream = 100000,
    downstream = 0,
    only_tss = FALSE,
    parallel = FALSE,
    tf_cor = 0.1,
    peak_cor = 0.,
    method = 'glm',
    alpha = 0.5,
    family = gaussian,
    interaction_term = ':',
    verbose = TRUE
){
    params <- NetworkParams(object)
    motif2tf <- NetworkTFs(object)
    regions <- NetworkRegions(object)
    features <- NetworkFeatures(object)
    gene_annot <- Signac::Annotation(object)
    gene_annot <- gene_annot[gene_annot$gene_name%in%features, ]

    log_message('Selecting candidate regulatory regions near genes', verbose=verbose)
    peak_data <- t(Seurat::GetAssayData(object, assay=params$peak_assay)[regions@peaks, ])
    colnames(peak_data) <- rownames(regions@motifs@data)
    peaks2motif <- regions@motifs@data

    # Find candidate regions near gene bodies
    peaks_near_gene <- find_peaks_near_genes(
        peaks = regions@ranges,
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
            model_frml, data=model_mat,
            family=family, method=method, alpha=alpha
        ), silent=TRUE)
        if (any(class(fit)=='try-error')){
            log_message('Warning: Fitting model failed for ', g, verbose=verbose==2)
            return()
        } else {
            fit$gof$nvariables <- nfeats
            return(fit)
        }
    }, verbose=verbose, parallel=parallel)
    model_fits <- model_fits[!map_lgl(model_fits, is.null)]
    coefs <- map_dfr(model_fits, function(x) x$coefs, .id='target')
    gof <- map_dfr(model_fits, function(x) x$gof, .id='target')

    network_obj <- new(
        Class = 'Network',
        coefs = coefs,
        fit = gof
    )
    object@grn@network <- network_obj
    object@grn@params[['method']] <- method
    object@grn@params[['family']] <- family
    object@grn@params[['dist']] <- c(upstream, downstream)
    object@grn@params[['only_tss']] <- only_tss
    object@grn@params[['interaction']] <- interaction_term
    object@grn@params[['tf_cor']] <- tf_cor
    object@grn@params[['peak_cor']] <- peak_cor

    return(object)

}



