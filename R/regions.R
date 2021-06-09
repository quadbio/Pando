#' @import dplyr
#' @importFrom Signac StringToGRanges
#' @importFrom Seurat GetAssay VariableFeatures
#' @importFrom S4Vectors subjectHits
#' @importFrom IRanges findOverlaps
NULL


#' Initiate the \code{RegulatoryNetwork} object.
#'
#' @param regions Candidate regions to consider for binding site inference.
#' If \code{NULL}, all peaks regions are considered.
#' @param genes A character vector with the genes to consider for GRN inference.
#' @param peak_assay A character vector indicating the name of the chromatin
#' accessibility assay in the \code{Seurat} object.
#' @param rna_assay A character vector indicating the name of the gene expression
#' assay in the \code{Seurat} object.
#' @param exclude_exons Logical. Whether to consider exons for binding site inference.
#' @param grn_name Character vector indicating the name of the GRN.
#'
#' @return A SeuratPlus object containing a RegulatoryNetwork object.
#'
#' @rdname initiate_grn
#' @export
#' @method initiate_grn Seurat
initiate_grn.Seurat <- function(
    object,
    regions = NULL,
    genes = VariableFeatures(object, assay=rna_assay),
    peak_assay = 'peaks',
    rna_assay = 'RNA',
    exclude_exons = TRUE
){
    peak_ranges <- StringToGRanges(rownames(GetAssay(object, assay=peak_assay)))

    if (!is.null(regions)){
        cand_ranges <- IRanges::intersect(regions, peak_ranges)
    } else {
        cand_ranges <- peak_ranges
    }

    peak_overlaps <- findOverlaps(cand_ranges, peak_ranges)
    peak_matches <- subjectHits(peak_overlaps)

    regions_obj <- new(
        Class = 'Regions',
        ranges = cand_ranges,
        peaks = peak_matches,
        motifs = NULL
    )

    gene_annot <- Signac::Annotation(object[[peak_assay]])
    genes_use <- intersect(gene_annot$gene_name, genes) %>%
        intersect(rownames(GetAssay(object, rna_assay)))
    genes <- list(
        genes = genes_use,
        tfs = NULL
    )

    params <- list(
        peak_assay = peak_assay,
        rna_assay = rna_assay,
        exclude_exons = exclude_exons
    )

    grn_obj <- new(
        Class = 'RegulatoryNetwork',
        regions = regions_obj,
        genes = genes,
        params = params
    )

    object <- as(object, 'SeuratPlus')
    object@grn <- grn_obj
    return(object)
}


#' Scan for motifs in candidate regions.
#'
#' @import dplyr
#'
#' @param pfm A \code{PFMatrixList} object with position weight matrices.
#' @param genome A \code{BSgenome} object with the genome of interest.
#' @param verbose Display messages.
#'
#' @return A SeuratPlus object with updated motif info.
#'
#' @rdname find_motifs
#' @export
#' @method find_motifs SeuratPlus
find_motifs.SeuratPlus <- function(
    object,
    pfm,
    genome,
    motif_tfs = NULL,
    verbose = TRUE
){
    params <- NetworkParams(object)

    # Add TF info for motifs
    log_message('Adding TF info', verbose=verbose)
    if (!is.null(motif_tfs)){
        motif2tf <- motif_tfs
    } else {
        utils::data(motif2tf, envir = environment())
    }

    motif2tf <- motif2tf %>% select(1,2) %>%
        distinct() %>% mutate(val=1) %>%
        tidyr::pivot_wider(names_from = 'tf', values_from=val, values_fill=0) %>%
        column_to_rownames('motif') %>% as.matrix() %>% Matrix::Matrix(sparse=TRUE)
    tfs_use <- intersect(rownames(GetAssay(object, params$rna_assay)), colnames(motif2tf))
    object@grn@genes$tfs <- motif2tf[, tfs_use]

    # Find motif positions with Signac/motifmatchr
    cand_ranges <- object@grn@regions@ranges
    motif_pos <- Signac::AddMotifs(
        object = cand_ranges,
        genome = genome,
        pfm = pfm,
        verbose= verbose
    )
    object@grn@regions@motifs <- motif_pos

    return(object)
}


