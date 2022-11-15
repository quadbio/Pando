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
#' @param methos Character specifying the method to link peak overlapping motif regions to nearby genes.
#' On of 'Signac' or 'GREAT'.
#' @param upstream Integer defining the distance upstream of the gene/TSS to consider.
#' @param downstream Integer defining the distance downstream of the gene/TSS to consider.
#' @param extend Integer defining the distance from the upstream and downstream of the
#' basal regulatory region. Used only when method is 'GREAT'.
#' @param sep Vector of separators to use for genomic string.
#' First element is used to separate chromosome and coordinates,
#' second separator is used to separate start and end coordinates.
#' @param only_tss Logical. Measure distance from the TSS (\code{TRUE})
#' or from the entire gene body (\code{FALSE}).
#' @param verbose Logical. Display messages
#'
#' @return A sparse binary Matrix with gene/peak matches.
#'
#' @export
find_peaks_near_genes <- function(
    peaks,
    genes,
    sep = c('-', '-'),
    method = c('Signac', 'GREAT'),
    upstream = 100000,
    downstream = 0,
    extend = 1000000,
    only_tss = FALSE,
    verbose = TRUE
){
    # Match arg
    method <- match.arg(method)

    if (method=='Signac'){

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

    } else if (method=='GREAT'){

        # Read gene annotation (Ensembl v93, GRCh38)
        utils::data(EnsDb.Hsapiens.v93.annot.UCSC.hg38, envir=environment())
        gene_annot_use <- EnsDb.Hsapiens.v93.annot.UCSC.hg38[
            which(EnsDb.Hsapiens.v93.annot.UCSC.hg38$gene_name %in% genes$gene_name),
        ]
        gene_annot_tss <- select(as_tibble(gene_annot_use), seqnames, 'start'=tss, 'end'=tss, strand)

        # Create GRanges object storing the TSS information
        tss <- GRanges(gene_annot_use)

        # Define basal regulatory region (promoter region)
        # as 5 kb upstream + 1 kb downstream of the TSS
        basal_reg <- suppressWarnings(
            expr = Signac::Extend(
                tss, upstream = upstream, downstream = downstream
            )
        )

        # Step 1 - get peaks overlap with basal regulatory region
        basal_overlaps <- suppressWarnings(IRanges::findOverlaps(
            query = peaks,
            subject = basal_reg,
            type = 'any',
            select = 'all',
            minoverlap = 2
        ))

        peak_all <- Signac::GRangesToString(grange = peaks, sep = sep)
        basal_peak_mapped_idx <- queryHits(basal_overlaps)
        basal_mapped_peaks <- unique(peak_all[basal_peak_mapped_idx])
        n1 <- length(basal_mapped_peaks)

        # Step 2: for the peaks not overlapped with basal regulatory regions,
        # check whether they located within gene body of any genes
        peak_unmapped_idx <- setdiff(seq(length(peak_all)), basal_peak_mapped_idx)
        peak_unmapped <- peak_all[peak_unmapped_idx]
        peak_unmapped_region <- Signac::StringToGRanges(peak_unmapped)

        # Create GRanges object storing annotated gene boundary
        gene_bound <- GRanges(gene_annot_use)
        body_overlaps <- IRanges::findOverlaps(
            query = peak_unmapped_region,
            subject = gene_bound,
            type = 'any',
            select = 'all',
            minoverlap = 2
        )
        body_peak_mapped_idx <- peak_unmapped_idx[queryHits(body_overlaps)]
        body_mapped_peaks <- unique(peak_all[body_peak_mapped_idx])
        n2 <- length(body_mapped_peaks)
        peak_mapped_idx <- c(basal_peak_mapped_idx, body_peak_mapped_idx)

        # Step 3: for the peaks not overlapped with basal regulatory regions of any genes,
        # check whether they overlap with extended regulatory region. i.e. +/- 1MB of basal regulatory region
        peak_unmapped_idx <- setdiff(seq(length(peak_all)), peak_mapped_idx)
        peak_unmapped <- peak_all[peak_unmapped_idx]
        peak_unmapped_region <- Signac::StringToGRanges(peak_unmapped)
        extend_reg <- suppressWarnings(
            expr = Signac::Extend(
                basal_reg, upstream = extend, downstream = extend
            )
        )

        # Get overlap between unmapped_peak_region and extended regulatory region
        extended_overlaps <- suppressWarnings(IRanges::findOverlaps(
            query = peak_unmapped_region,
            subject = extend_reg,
            type = 'any',
            select = 'all',
            minoverlap = 2
        ))
        extended_peak_mapped_idx <- peak_unmapped_idx[queryHits(extended_overlaps)]
        extended_mapped_peaks <- unique(peak_all[extended_peak_mapped_idx])
        n3 <- length(extended_mapped_peaks)

        hit_matrix <- Matrix::sparseMatrix(
            i = c(basal_peak_mapped_idx,
                  body_peak_mapped_idx,
                  extended_peak_mapped_idx),
            j = c(subjectHits(basal_overlaps),
                  subjectHits(body_overlaps),
                  subjectHits(extended_overlaps)),
            x = 1,
            dims = c(length(peaks), length(basal_reg))
        )
        rownames(hit_matrix) <- peak_all
        colnames(hit_matrix) <- c(basal_reg$gene_name)
    }
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

#' Copy of the aggregate.Matrix function from the Matrix.utils package,
#' since this is off CRAN and does not seem to be maintained anymore
#' @keyword internal
#'
fast_aggregate <- function(
    x,
    groupings = NULL,
    form = NULL,
    fun = 'sum',
    ...
){
    if (!is(x,'Matrix')){
        x <- Matrix(as.matrix(x), sparse=TRUE)
    }
    if (fun=='count'){
        x <- x!=0
    }
    groupings2 <- groupings
    if (!is(groupings2, 'data.frame')){
        groupings2 <- as.data.frame(groupings2)
    }
    groupings2 <- data.frame(lapply(groupings2, as.factor))
    groupings2 <- data.frame(interaction(groupings2, sep='_'))
    colnames(groupings2) <- 'A'
    if (is.null(form)){
        form <- as.formula('~0+.')
    }
    form <- as.formula(form)
    mapping <- dMcast(groupings2, form)
    colnames(mapping) <- substring(colnames(mapping), 2)
    result <- t(mapping) %*% x
    if (fun=='mean'){
        result@x <- result@x/(fast_aggregate(x, groupings2, fun='count'))@x
    }
    attr(result,'crosswalk') <- grr::extract(groupings, match(rownames(result), groupings2$A))
    return(result)
}

#' Copy of the dMcast function from the Matrix.utils package,
#' since this is off CRAN and does not seem to be maintained anymore
#' @keyword internal
#'
dMcast <- function(
    data,
    formula,
    fun.aggregate = 'sum',
    value.var = NULL,
    as.factors = FALSE,
    factor.nas = TRUE,
    drop.unused.levels = TRUE
){
    values <- 1
    if (!is.null(value.var)){
        values <- data[,value.var]
    }
    alltms <- terms(formula, data=data)
    response <- rownames(attr(alltms, 'factors'))[attr(alltms, 'response')]
    tm <- attr(alltms, "term.labels")
    interactionsIndex <- grep(':', tm)
    interactions <- tm[interactionsIndex]
    simple <- setdiff(tm, interactions)
    i2 <- strsplit(interactions,':')
    newterms <- unlist(lapply(i2, function (x) paste("paste(", paste(x, collapse=','), ",", "sep='_'",")")))
    newterms <- c(simple, newterms)
    newformula <- as.formula(paste('~0+', paste(newterms, collapse='+')))
    allvars <- all.vars(alltms)
    data <- data[, c(allvars), drop=FALSE]
    if (as.factors)
        data <- data.frame(lapply(data, as.factor))
    characters <- unlist(lapply(data, is.character))
    data[,characters] <- lapply(data[, characters,drop=FALSE], as.factor)
    factors <- unlist(lapply(data, is.factor))
    # Prevents errors with 1 or fewer distinct levels
    data[,factors] <- lapply(data[,factors,drop=FALSE],function (x)
    {
        if (factor.nas){
            if (any(is.na(x))){
                levels(x) <- c(levels(x),'NA')
                x[is.na(x)] <- 'NA'
            }
        }
        if (drop.unused.levels){
            if (nlevels(x)!=length(na.omit(unique(x)))){
                x <- factor(as.character(x))
            }
        }
        y <- contrasts(x, contrasts=FALSE, sparse=TRUE)
        attr(x, 'contrasts') <- y
        return(x)
    })
    # Allows NAs to pass
    attr(data,'na.action') <- na.pass
    result <- Matrix::sparse.model.matrix(newformula, data,drop.unused.levels = FALSE, row.names=FALSE)
    brokenNames <- grep('paste(', colnames(result), fixed = TRUE)
    colnames(result)[brokenNames] <- lapply(colnames(result)[brokenNames], function (x) {
        x <- gsub('paste(', replacement='', x=x, fixed = TRUE)
        x <- gsub(pattern=', ', replacement='_', x=x, fixed=TRUE)
        x <- gsub(pattern='_sep = \"_\")', replacement='', x=x, fixed=TRUE)
        return(x)
    })

    result <- result*values
    if(isTRUE(response>0))
    {
        responses=all.vars(terms(as.formula(paste(response,'~0'))))
        result <- fast_aggregate(result, data[, responses,drop=FALSE], fun=fun.aggregate)
    }
    return(result)
}


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
    if (length(groups) == nrow(x) & 'character'%in%class(fun)){
        if (fun%in%c('count', 'sum')){
            agg_mat <- fast_aggregate(x=x, groupings=groups, fun=fun)
            return(agg_mat)
        }

        if (fun=='mean'){
            group_counts <- as.numeric(table(groups))
            agg_mat <- fast_aggregate(x=x, groupings=groups, fun='sum')
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
        agg_mat <- Matrix::Matrix(agg_mat, sparse=TRUE)
    } else if (length(groups) <= 1){
        agg_mat <- fun(x)
        agg_mat <- Matrix::Matrix(agg_mat, sparse=TRUE)
        colnames(agg_mat) <- groups
        rownames(agg_mat) <- colnames(x)
    } else {
        stop('Length of groups must be either nrow(x) or 1.')
    }
    return(Matrix::t(agg_mat))
}


#' Aggregate Seurat assay over groups
#'
#' @param group_name A character vector indicating the metadata column to aggregate over.
#' @param fun The summary function to be applied to each group.
#' @param assay The assay to summarize.
#' @param slot The slot to summarize.
#'
#' @return A Seurat object.
#'
#' @export
aggregate_assay <- function(
    object,
    group_name,
    fun = 'mean',
    assay = 'RNA',
    slot = 'data'
){
    ass_mat <- Matrix::t(Seurat::GetAssayData(object, assay=assay, slot=slot))
    groups <- as.character(object@meta.data[[group_name]])
    agg_mat <- aggregate_matrix(ass_mat, groups=groups, fun=fun)
    if (is.null(object@assays[[assay]]@misc$summary)){
        object@assays[[assay]]@misc$summary <- list()
    }
    object@assays[[assay]]@misc$summary[[group_name]] <- agg_mat
    return(object)
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
    } else if (method == 'spearman' & require(presto, quietly=TRUE)){
        xr <- Matrix(presto::rank_matrix(t(x))$X_ranked, sparse=TRUE)
        rownames(xr) <- rownames(x)
        colnames(xr) <- colnames(x)
        if (!is.null(y)){
            yr <- Matrix(presto::rank_matrix(t(y))$X_ranked, sparse=TRUE)
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

