#' @importFrom Signac StringToGRanges
#' @importFrom Seurat GetAssay VariableFeatures
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom IRanges findOverlaps pintersect
NULL


#' Initiate the \code{RegulatoryNetwork} object.
#'
#' @param regions Candidate regions to consider for binding site inference.
#' If \code{NULL}, all peaks regions are considered.
#' @param peak_assay A character vector indicating the name of the chromatin
#' accessibility assay in the \code{Seurat} object.
#' @param rna_assay A character vector indicating the name of the gene expression
#' assay in the \code{Seurat} object.
#' @param exclude_exons Logical. Whether to consider exons for binding site inference.
#'
#' @return A GRNData object containing a RegulatoryNetwork object.
#'
#' @rdname initiate_grn
#' @export
#' @method initiate_grn Seurat
initiate_grn.Seurat <- function(
    object,
    regions = NULL,
    peak_assay = 'peaks',
    rna_assay = 'RNA',
    exclude_exons = TRUE
){
    gene_annot <- Signac::Annotation(object[[peak_assay]])
    # Through error if NULL
    if (is.null(gene_annot)){
        stop('Please provide a gene annotation for the ChromatinAssay.')
    }
    peak_ranges <- Signac::StringToGRanges(
        rownames(Seurat::GetAssay(object, assay=peak_assay))
    )

    # Find candidate ranges by intersecting the supplied regions with peaks
    # Per default take all peaks as candidate regions
    if (!is.null(regions)){
        cand_olaps <- findOverlaps(regions, peak_ranges)
        cand_ranges <- pintersect(
            peak_ranges[subjectHits(cand_olaps)],
            regions[queryHits(cand_olaps)]
        )
    } else {
        cand_ranges <- peak_ranges
    }

    # Exclude exons because they are usually conserved
    if (exclude_exons){
        exon_ranges <- gene_annot[gene_annot$type=='exon', ]
        names(exon_ranges@ranges) <- NULL
        exon_ranges <- IRanges::intersect(exon_ranges, exon_ranges)
        exon_ranges <- GenomicRanges::GRanges(
            seqnames = exon_ranges@seqnames,
            ranges = exon_ranges@ranges
        )
        cand_ranges <- GenomicRanges::subtract(
            cand_ranges, exon_ranges, ignore.strand=TRUE
        ) %>% unlist()
    }

    # Match candidate ranges to peaks
    peak_overlaps <- findOverlaps(cand_ranges, peak_ranges)
    peak_matches <- subjectHits(peak_overlaps)

    regions_obj <- new(
        Class = 'Regions',
        ranges = cand_ranges,
        peaks = peak_matches,
        motifs = NULL
    )

    params <- list(
        peak_assay = peak_assay,
        rna_assay = rna_assay,
        exclude_exons = exclude_exons
    )

    # Add the modified peak matches to the Seurat object
    cand_regions_df <- as.data.frame(cand_ranges) %>% 
                      unite('region', c('seqnames', 'start', 'end'), sep='-', remove=FALSE) %>% 
                      column_to_rownames('region')
    cand_regions <- rownames(cand_regions_df)

    cand_assay <- object[[peak_assay]]
    cand_assay@counts <- cand_assay@counts[peak_matches, ]
    rownames(cand_assay@counts) <- cand_regions
    cand_assay@data <- cand_assay@data[peak_matches, ]
    rownames(cand_assay@data) <- cand_regions
    cand_assay@ranges <- cand_ranges
    cand_assay@var.features <- list()
    cand_assay@meta.features <- cand_regions_df

    object[[paste0(peak_assay, "_cand")]] <- cand_assay

    grn_obj <- new(
        Class = 'RegulatoryNetwork',
        regions = regions_obj,
        params = params
    )

    object <- new(
        Class = 'GRNData',
        grn = grn_obj,
        data = object
    )

    return(object)
}

#' Initiate the \code{RegulatoryNetwork} object.
#'
#' @param object An object.
#'
#' @rdname initiate_grn
#' @export
#' @method initiate_grn GRNData
initiate_grn.GRNData <- function(
        object,
        regions = NULL,
        peak_assay = 'peaks',
        rna_assay = 'RNA',
        exclude_exons = TRUE
){
    initiate_grn(object@data)
}

#' Scan for motifs in candidate regions.
#'
#' @importFrom dplyr distinct mutate select
#'
#' @param pfm A \code{PFMatrixList} object with position weight matrices.
#' @param genome A \code{BSgenome} object with the genome of interest.
#' @param motif_tfs A data frame matching motifs with TFs. The first column is assumed
#' to be the name of the motif, the second the name of the TF.
#' @param verbose Display messages.
#'
#' @return A GRNData object with updated motif info.
#'
#' @rdname find_motifs
#' @export
#' @method find_motifs GRNData
find_motifs.GRNData <- function(
    object,
    pfm,
    genome,
    motif_tfs = NULL,
    verbose = TRUE
){
    params <- Params(object)

    # Add TF info for motifs
    log_message('Adding TF info', verbose=verbose)
    if (!is.null(motif_tfs)){
        motif2tf <- motif_tfs
    } else {
        utils::data(motif2tf, envir = environment())
    }

    # Spread dataframe to sparse matrix
    motif2tf <- motif2tf %>% select('motif'=1,'tf'=2) %>%
        distinct() %>% mutate(val=1) %>%
        tidyr::pivot_wider(names_from = 'tf', values_from=val, values_fill=0) %>%
        tibble::column_to_rownames('motif') %>%
        as.matrix() %>% Matrix::Matrix(sparse=TRUE)
    tfs_use <- intersect(
        rownames(GetAssay(object, params$rna_assay)),
        colnames(motif2tf)
    )
    if (length(tfs_use)==0){
        stop('None of the provided TFs were found in the dataset. Consider providing a custom motif-to-TF map as `motif_tfs`')
    }
    object@grn@regions@tfs <- motif2tf[, tfs_use]

    # Find motif positions with Signac/motifmatchr
    cand_ranges <- object@grn@regions@ranges
    motif_pos <- Signac::AddMotifs(
        object = cand_ranges,
        genome = genome,
        pfm = pfm,
        verbose= verbose
    )
    # Add motif positions to grn object
    object@grn@regions@motifs <- motif_pos
    # Also add motif positions to the data object for default assay
    object@data[[paste0(params$peak_assay, "_cand")]]@motifs <- motif_pos

    return(object)
}
