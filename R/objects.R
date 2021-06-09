#' @importFrom methods setClass
#' @importClassesFrom Signac Motif
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom SeuratObject Seurat
NULL


#' The Network class
#'
#' The Network object stores the inferred network itself, information about the fitting
#' process as well as graph representations of the network.
#'
#' @slot fit A dataframe with goodness of fit measures.
#' @slot coefs A dataframe with the fitted coefficients.
#' @slot graph A graphical representation of the inferred network.
#'
#' @name Network-class
#' @rdname Network-class
#' @exportClass Network
Network <- setClass(
    Class = 'Network',
    slots = list(
        fit = 'data.frame',
        coefs = 'data.frame',
        graph = 'ANY'
    )
)


#' The Regions class
#'
#' The Regions object stores the genomic regions that are considered by the model.
#' It stores their genomic positions, how they map to the peaks in the Seurat object
#' and motif matches.
#'
#' @slot motifs A \code{Motifs} object with matches of TF motifs.
#' @slot ranges A \code{GenomicRanges} object.
#' @slot peaks A numeric vector with peak indices for each region.
#'
#' @name Regions-class
#' @rdname Regions-class
#' @exportClass Regions
Regions <- setClass(
    Class = 'Regions',
    slots = list(
        motifs = 'ANY',
        ranges = 'GRanges',
        peaks = 'numeric'
    )
)


#' The RegularotyNetwork class
#'
#' The RegularotyNetwork object is the core data structure in Pando.
#' It stores all data necessary for network inference and analysis
#' that is not provided by Seurat.
#'
#' @slot genes A named list containing the transcription factors and
#' target genes included in the network.
#' @slot regions A \code{\link{Regions}} object containing information about
#' the genomic regions included in the network.
#' @slot network A \code{\link{Network}} object containing the inferred regulatory
#' network and information about the model fit.
#' @slot params A list storing parameters for GRN inference.
#'
#' @name RegulatoryNetwork-class
#' @rdname RegulatoryNetwork-class
#' @exportClass RegulatoryNetwork
RegulatoryNetwork <- setClass(
    Class = 'RegulatoryNetwork',
    slots = list(
        genes = 'list',
        regions = 'Regions',
        network = 'Network',
        params = 'list'
    )
)


#' The SeuratPlus class
#'
#' The SeuratPlus object is an extended \code{Seurat} object
#' for the storage and analysis of Regulatory network data.
#'
#' @slot grn A named list containing \code{RegulatoryNetwork} objects with inferred networks.
#'
#' @name SeuratPlus-class
#' @rdname SeuratPlus-class
#' @exportClass SeuratPlus
#' @concept assay
SeuratPlus <- setClass(
    Class = 'SeuratPlus',
    contains = 'Seurat',
    slots = list(
        'grn' = 'RegulatoryNetwork'
    )
)


#' Get network features
#' @rdname NetworkFeatures
#' @method NetworkFeatures SeuratPlus
#' @export
NetworkFeatures.SeuratPlus <- function(object){
    return(object@grn@genes$genes)
}

#' @rdname NetworkFeatures
#' @method NetworkFeatures RegulatoryNetwork
#' @export
NetworkFeatures.RegulatoryNetwork <- function(object){
    return(object@genes$genes)
}


#' Get network TFs
#' @rdname NetworkTFs
#' @method NetworkTFs SeuratPlus
#' @export
NetworkTFs.SeuratPlus <- function(object){
    return(object@grn@genes$tfs)
}

#' @rdname NetworkTFs
#' @method NetworkTFs RegulatoryNetwork
#' @export
NetworkTFs.RegulatoryNetwork <- function(object){
    return(object@genes$tfs)
}


#' Get network regions
#' @rdname NetworkRegions
#' @method NetworkRegions SeuratPlus
#' @export
NetworkRegions.SeuratPlus <- function(object){
    return(object@grn@regions)
}

#' @rdname NetworkRegions
#' @method NetworkRegions RegulatoryNetwork
#' @export
NetworkRegions.RegulatoryNetwork <- function(object){
    return(object@regions)
}


#' Get network data
#' @rdname GetNetworkData
#' @method GetNetworkData SeuratPlus
#' @export
GetNetworkData.SeuratPlus <- function(object){
    return(object@grn)
}


#' Get fitted coefficients
#' @rdname coef
#' @method coef SeuratPlus
#' @export
coef.SeuratPlus <- function(object){
    return(object@grn@network@coefs)
}

#' @rdname coef
#' @method coef RegulatoryNetwork
#' @export
coef.RegulatoryNetwork <- function(object){
    return(object@network@coefs)
}


#' Get goodness-of-fit info
#' @rdname gof
#' @method gof SeuratPlus
#' @export
gof.SeuratPlus <- function(object){
    return(object@grn@network@gof)
}

#' @rdname gof
#' @method gof RegulatoryNetwork
#' @export
gof.RegulatoryNetwork <- function(object){
    return(object@network@gof)
}


#' Get network parameters
#' @rdname NetworkParams
#' @method NetworkParams SeuratPlus
#' @export
NetworkParams.SeuratPlus <- function(object){
    return(object@grn@params)
}

#' @rdname NetworkParams
#' @method NetworkParams RegulatoryNetwork
#' @export
NetworkParams.RegulatoryNetwork <- function(object){
    return(object@params)
}


#' Print RegulatoryNetwork objects
#'
#' @rdname print
#' @export
#' @method print RegulatoryNetwork
print.RegulatoryNetwork <- function(object){
    n_genes <- length(NetworkFeatures(object))
    n_tfs <- ncol(NetworkTFs(object))
    if (is.null(n_tfs)){
        tf_string <- '\nCandidate regions have not been scanned for motifs'
    } else {
        tf_string <- paste0('\nand ', n_tfs, ' transcription factors')
    }
    n_conn <- coef(object)
    if (all(n_conn==0)){
        conn_string <- '\nNetwork has not been inferred'
    } else {
        conn_string <- NULL
    }
    cat(paste0(
        'An object of class RegulatoryNetwork\n', 'based on ', n_genes, ' target genes',
        tf_string, conn_string
    ))
}

setMethod('show', 'RegulatoryNetwork', function(object) print(object))






