#' Print diagnostic message
#'
#' @param ... Text to print
#' @param verbose Display messages
log_message <- function(..., verbose=T){
    if (verbose){
        message(paste0(...))
    }
}



#' Find peaks or regions near gene body or TSS
#'
#' @param peaks A \code{GRanges} object with peak regions.
#' @param gene A \code{GRanges} object with gene coordinates.
#' @param upstream Integer defining the distance upstream of the gene to consider.
#' @param downstream Integer defining the distance downstream of the gene to consider.
#' @param sep Vector of separators to use for genomic string.
#' First element is used to separate chromosome and coordinates,
#' second separator is used to separate start and end coordinates.
#' @param only_tss Logical. Measure distance from the TSS (\code{TRUE}) or from the entire gene body (\code{FALSE}).
#' @param verbose Logical. Display messages
#'
#' @return A sparse binary Matrix with gene/peak matches.
find_peaks_near_genes <- function(
    peaks,
    genes,
    upstream = 100000,
    downstream = 0,
    sep = c('-', '-'),
    only_tss = FALSE
){
    if (only_tss){
        genes <- IRanges::resize(x = genes, width = 1, fix = 'start')
    }
    genes_extended <- suppressWarnings(
        expr = Signac::Extend(
            genes, upstream = upstream, downstream = downstream
        )
    )
    overlaps <- IRanges::findOverlaps(
        query = peaks,
        subject = genes_extended,
        type = 'any',
        select = 'all'
    )
    hit_matrix <- Matrix::sparseMatrix(
        i = S4Vectors::queryHits(overlaps),
        j = S4Vectors::subjectHits(overlaps),
        x = 1,
        dims = c(length(peaks), length(genes_extended))
    )
    rownames(hit_matrix) <- Signac::GRangesToString(grange = peaks, sep = sep)
    colnames(hit_matrix) <- genes_extended$gene_name
    return(hit_matrix)
}



#' @import sparseMatrixStats
summary_fun <- list(
    'mean' = sparseMatrixStats::colMeans2,
    'median' = sparseMatrixStats::colMedians,
    'max' = sparseMatrixStats::colMaxs,
    'min' = sparseMatrixStats::colMins,
    'count' = sparseMatrixStats::colCounts,
    'any' = sparseMatrixStats::colAnys,
    'all' = sparseMatrixStats::colAlls,
    'sd' = sparseMatrixStats::colSds,
    'mad' = sparseMatrixStats::colMads
)



#' Aggregate matrix over groups
#'
#' @import sparseMatrixStats
#'
#' @param groups A character vector with the groups to aggregate over.
#' @param fun The summary function to be applied to each group.
#'
#' @return A summary matrix.
#'
#' @export
aggregate_matrix <- function(
    x,
    groups = NULL,
    fun = 'mean'
){

    if (length(groups) == nrow(x)){
        if (fun%in%c('count', 'sum')){
            agg_mat <- Matrix.utils::aggregate.Matrix(x=x, groupings=groups, fun=fun)
            return(agg_mat)
        }

        if (fun=='mean'){
            group_counts <- as.numeric(table(groups))
            agg_mat <- Matrix.utils::aggregate.Matrix(x=x, groupings=groups, fun='sum')
            agg_mat <- agg_mat / group_counts
            return(agg_mat)
        }
    }

    if ('character'%in%class(fun)){
        fun <- summary_fun[[fun]]
    }

    if (length(groups) == nrow(x)){
        agg_mat <- sapply(levels(factor(groups)), function(g){
            chunk <- x[which(groups==g), ]
            if (is.null(dim(chunk))){
                return(chunk)
            } else {
                return(fun(chunk))
            }
        })
        agg_mat <- Matrix::Matrix(agg_mat, sparse=T)
    } else if (length(groups) <= 1){
        agg_mat <- fun(x)
        agg_mat <- Matrix::Matrix(agg_mat, sparse=T)
        colnames(agg_mat) <- groups
        rownames(agg_mat) <- colnames(x)
    } else {
        stop('Length of groups must be either nrow(x) or 1.')
    }
    return(Matrix::t(agg_mat))
}


#' Apply function over a List or Vector (in parallel)
#'
#' @param x A vector or list to apply over.
#' @param fun The function to be applied to each element.
#' @param parallel Logical. Whether to use \code{foreach} to parallelize.
#' @param verbose Logical. Whether to print progress bar.
#' Only works in sequential mode.
#'
#' @return A list.
#'
#' @export
map_par <- function(x, fun, parallel=FALSE, verbose=TRUE){
    if (!parallel & (verbose==1)){
        return(pbapply::pblapply(X=x, FUN=fun))
    }
    if (!parallel & (verbose!=1)){
        return(base::lapply(X=x, FUN=fun))
    }
    if (parallel){
        outlist <- foreach::foreach(i=1:length(x)) %dopar% {fun(x[[i]])}
        names(outlist) <- names(x)
        return(outlist)
    }
}



#' Fast correlation and covariance calcualtion for sparse matrices
#'
#' @import Matrix
#'
#' @param x Sparse matrix or character vector.
#' @param y Sparse matrix or character vector.
#'
#' @return A list containing a covariance and correlation matrix.
sparse_covcor <- function(x, y=NULL) {
    if (!is(x, "dgCMatrix")) stop("x should be a dgCMatrix")
    if (is.null(y)) {
        n <- nrow(x)
        muX <- colMeans(x)
        covmat <- (as.matrix(crossprod(x)) - n * tcrossprod(muX)) / (n-1)
        sdvec <- sqrt(diag(covmat))
        cormat <- covmat / tcrossprod(sdvec)
        return(list(cov=covmat, cor=cormat))
    } else {
        if (!is(y, "dgCMatrix")) stop("y should be a dgCMatrix")
        if (nrow(x) != nrow(y)) stop("x and y should have the same number of rows")
        n <- nrow(x)
        muY <- colMeans(y)
        muX <- colMeans(x)
        covmat <- (as.matrix(crossprod(x, y)) - n * tcrossprod(muX, muY)) / (n-1)
        sdvecX <- sqrt((colSums(x^2) - n*muX^2) / (n-1))
        sdvecY <- sqrt((colSums(y^2) - n*muY^2) / (n-1))
        cormat <- covmat / tcrossprod(sdvecX, sdvecY)
        return(list(cov=covmat, cor=cormat))
    }
}


#' Safe correlation function which returns a sparse matrix without missing values
#'
#' @import Matrix
#'
#' @param x Sparse matrix or character vector.
#' @param y Sparse matrix or character vector.
#' @param method Method to use for calculating the correlation coefficient.
#' @param allow_neg Logical. Whether to allow negative values or set them to 0.
#' @param ... Other arguments passed to the correlation function.
#'
#' @return A correlation matrix.
#'
#' @export
sparse_cor <- function(
    x,
    y = NULL,
    method = 'pearson',
    allow_neg = TRUE,
    remove_na = TRUE,
    remove_inf = TRUE,
    ...
){
    if (method == 'pearson'){
        x <- Matrix(x, sparse=TRUE)
        if (!is.null(y)){
            y <- Matrix(y, sparse=TRUE)
        }
        corr_mat <- sparse_covcor(x, y)$cor
    } else if (method == 'spearman' & require(presto, quietly=T)){
        xr <- Matrix(presto::rank_matrix(t(x))$X_ranked, sparse=T)
        rownames(xr) <- rownames(x)
        colnames(xr) <- colnames(x)
        if (!is.null(y)){
            yr <- Matrix(presto::rank_matrix(t(y))$X_ranked, sparse=T)
            rownames(yr) <- rownames(y)
            colnames(yr) <- colnames(y)
            y <- yr
        }
        corr_mat <- sparse_covcor(xr, y)$cor
    } else {
        x <- as.matrix(x)
        if (!is.null(y)){
            y <- as.matrix(y)
        }
        corr_mat <- stats::cor(x, y, method=method, ...)
    }
    if (remove_na){
        corr_mat[is.na(corr_mat)] <- 0
    }
    if (remove_inf){
        corr_mat[is.infinite(corr_mat)] <- 1
    }
    corr_mat <- Matrix(corr_mat, sparse=TRUE)
    if (!allow_neg){
        corr_mat[corr_mat < 0] <- 0
    }
    return(corr_mat)
}
