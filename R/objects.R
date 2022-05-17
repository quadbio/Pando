#' @importFrom methods setClass
#' @importClassesFrom Signac Motif
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom SeuratObject Seurat
NULL


#' The Modules class
#'
#' The Modules object stores the TF modules extracted from the inferred network..
#'
#' @slot meta A dataframe with meta data about the modules.
#' @slot features A named list with a set of fetures (genes/peaks) for each module.
#' @slot params A named list with module selection parameters.
#'
#' @name Network-class
#' @rdname Network-class
#' @exportClass Network
Modules <- setClass(
    Class = 'Modules',
    slots = list(
        meta = 'data.frame',
        features = 'list',
        params = 'list'
    )
)


#' The Network class
#'
#' The Network object stores the inferred network itself, information about the fitting
#' process as well as graph representations of the network.
#'
#' @slot fit A dataframe with goodness of fit measures.
#' @slot coefs A dataframe with the fitted coefficients.
#' @slot modules A list TF modules.
#' @slot graphs Graphical representations of the inferred network.
#' @slot params A named list with GRN inference parameters.
#'
#' @name Network-class
#' @rdname Network-class
#' @exportClass Network
Network <- setClass(
    Class = 'Network',
    slots = list(
        features = 'character',
        fit = 'data.frame',
        coefs = 'data.frame',
        modules = 'Modules',
        graphs = 'list',
        params = 'list'
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
        tfs = 'ANY',
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
#' @slot active_network A string indicating the active network.
#'
#' @name RegulatoryNetwork-class
#' @rdname RegulatoryNetwork-class
#' @exportClass RegulatoryNetwork
RegulatoryNetwork <- setClass(
    Class = 'RegulatoryNetwork',
    slots = list(
        regions = 'Regions',
        networks = 'list',
        params = 'list',
        active_network = 'character'
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

#' Get network
#' @rdname GetNetwork
#' @method GetNetwork SeuratPlus
#' @export
GetNetwork.SeuratPlus <- function(object, network=DefaultNetwork(object)){
    return(GetNetwork(object@grn, network=network))
}

#' Get network
#' @rdname GetNetwork
#' @method GetNetwork RegulatoryNetwork
#' @export
GetNetwork.RegulatoryNetwork <- function(object, network=DefaultNetwork(object)){
    if (length(network)==0){
        stop(paste0('Network not found, try running `infer_network()` first.'))
    }
    if (!network%in%names(object@networks)){
        stop(paste0('The requested network "', network, '" does not exist.'))
    }
    return(object@networks[[network]])
}


#' Get network features
#' @rdname NetworkFeatures
#' @method NetworkFeatures SeuratPlus
#' @export
NetworkFeatures.SeuratPlus <- function(object, network=DefaultNetwork(object)){
    return(GetNetwork(object, network=network)@features)
}

#' @rdname NetworkFeatures
#' @method NetworkFeatures RegulatoryNetwork
#' @export
NetworkFeatures.RegulatoryNetwork <- function(object, network=DefaultNetwork(object)){
    return(GetNetwork(object, network=network)@features)
}


#' Get network TFs
#' @rdname NetworkTFs
#' @method NetworkTFs SeuratPlus
#' @export
NetworkTFs.SeuratPlus <- function(object){
    return(object@grn@regions@tfs)
}

#' @rdname NetworkTFs
#' @method NetworkTFs RegulatoryNetwork
#' @export
NetworkTFs.RegulatoryNetwork <- function(object){
    return(object@regions@tfs)
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
#' @rdname GetGRN
#' @method GetGRN SeuratPlus
#' @export
GetGRN.SeuratPlus <- function(object){
    return(object@grn)
}


#' Get TF modules
#' @rdname NetworkModules
#' @method NetworkModules SeuratPlus
#' @export
NetworkModules.SeuratPlus <- function(object, network=DefaultNetwork(object)){
    return(GetNetwork(object, network=network)@modules)
}

#' @rdname NetworkModules
#' @method NetworkModules RegulatoryNetwork
#' @export
NetworkModules.RegulatoryNetwork <- function(object, network=DefaultNetwork(object)){
    return(GetNetwork(object, network=network)@modules)
}

#' @rdname NetworkModules
#' @method NetworkModules Network
#' @export
NetworkModules.Network <- function(object){
    return(object@modules)
}


#' Get network parameters
#' @rdname NetworkParams
#' @method NetworkParams SeuratPlus
#' @export
NetworkParams.SeuratPlus <- function(object, network=DefaultNetwork(object)){
    return(GetNetwork(object, network=network)@params)
}

#' @rdname NetworkParams
#' @method NetworkParams RegulatoryNetwork
#' @export
NetworkParams.RegulatoryNetwork <- function(object, network=DefaultNetwork(object)){
    return(GetNetwork(object, network=network)@params)
}

#' @rdname NetworkParams
#' @method NetworkParams Network
#' @export
NetworkParams.Network <- function(object){
    return(object@params)
}


#' Get network parameters
#' @rdname NetworkGraph
#' @method NetworkGraph SeuratPlus
#' @export
NetworkGraph.SeuratPlus <- function(object, network=DefaultNetwork(object), graph='module_graph'){
    return(NetworkGraph(GetNetwork(object, network=network), graph=graph))
}

#' @rdname NetworkGraph
#' @method NetworkGraph RegulatoryNetwork
#' @export
NetworkGraph.RegulatoryNetwork <- function(object, network=DefaultNetwork(object), graph='module_graph'){
    return(NetworkGraph(GetNetwork(object, network=network), graph=graph))
}

#' @rdname NetworkGraph
#' @method NetworkGraph Network
#' @export
NetworkGraph.Network <- function(object, graph='module_graph'){
    if (!graph%in%names(object@graphs)){
        stop(paste0('The requested graph "', graph, '" does not exist. Try (re-)running `get_network_graph().`'))
    }
    return(object@graphs[[graph]])
}

#' Get active network
#' @rdname DefaultNetwork
#' @method DefaultNetwork SeuratPlus
#' @export
DefaultNetwork.SeuratPlus <- function(object){
    return(DefaultNetwork(GetGRN(object)))
}

#' @rdname DefaultNetwork
#' @method DefaultNetwork RegulatoryNetwork
#' @export
DefaultNetwork.RegulatoryNetwork <- function(object){
    return(object@active_network)
}


#' Get fitted coefficients
#' @rdname coef
#' @method coef SeuratPlus
#' @export
coef.SeuratPlus <- function(object, network=DefaultNetwork(object)){
    return(GetNetwork(object, network=network)@coefs)
}

#' @rdname coef
#' @method coef RegulatoryNetwork
#' @export
coef.RegulatoryNetwork <- function(object, network=DefaultNetwork(object)){
    return(GetNetwork(object, network=network)@coefs)
}

#' @rdname coef
#' @method coef Network
#' @export
coef.Network <- function(object){
    return(object@coefs)
}


#' Get goodness-of-fit info
#' @rdname gof
#' @method gof SeuratPlus
#' @export
gof.SeuratPlus <- function(object, network=DefaultNetwork(object)){
    return(GetNetwork(object, network=network)@fit)
}

#' @rdname gof
#' @method gof RegulatoryNetwork
#' @export
gof.RegulatoryNetwork <- function(object, network=DefaultNetwork(object)){
    return(GetNetwork(object, network=network)@fit)
}

#' @rdname gof
#' @method gof Network
#' @export
gof.Network <- function(object, network=DefaultNetwork(object)){
    return(object@fit)
}


#' Get GRN inference parameters
#' @rdname Params
#' @method Params SeuratPlus
#' @export
Params.SeuratPlus <- function(object){
    return(object@grn@params)
}

#' @rdname Params
#' @method Params RegulatoryNetwork
#' @export
Params.RegulatoryNetwork <- function(object){
    return(object@params)
}


#' Get summary of seurat assay
#' @rdname GetAssaySummary
#' @method GetAssaySummary Seurat
#' @export
GetAssaySummary.Seurat <- function(object, group_name, assay=NULL, verbose=TRUE){
    if (is.null(assay)){
        assay <- object@active.assay
    }
    smry <- Seurat::Misc(object[[assay]])$summary[[group_name]]
    if (is.null(smry)){
        log_message('Summary of "', group_name, '" does not yet exist.', verbose=verbose)
        log_message('Summarizing.', verbose=verbose)
        object <- aggregate_assay(object, assay=assay, group_name=group_name)
        smry <- GetAssaySummary(object, assay=assay, group_name=group_name, verbose=verbose)
    }
    return(smry)
}


#' Print RegulatoryNetwork objects
#'
#' @rdname print
#' @export
#' @method print RegulatoryNetwork
print.RegulatoryNetwork <- function(object){
    n_tfs <- ncol(NetworkTFs(object))
    if (is.null(n_tfs)){
        tf_string <- '\nCandidate regions have not been scanned for motifs'
    } else {
        tf_string <- paste0('based on ', n_tfs, ' transcription factors')
    }
    n_nets <- length(object@networks)
    net_names <- names(object@networks)
    if (n_nets==0){
        conn_string <- '\nNo network has been inferred'
    } else if (n_nets==1){
        conn_string <- paste0(n_nets, ' inferred network: ', net_names)
    } else {
        conn_string <- paste0(n_nets, ' inferred networks: ',
            paste(net_names, collapse=', '))
    }
    cat(paste0(
        'A RegulatoryNetwork object ', tf_string, '\n',
        conn_string
    ))
}

setMethod('show', 'RegulatoryNetwork', function(object) print(object))


#' Print Network objects
#'
#' @rdname print
#' @export
#' @method print Network
print.Network <- function(object){
    if (nrow(NetworkModules(object)@meta)==0){
        n_genes <- length(unique(coef(object)$target))
        n_tfs <- length(unique(coef(object)$tf))
    } else {
        n_genes <- length(unique(NetworkModules(object)@meta$target))
        n_tfs <- length(unique(NetworkModules(object)@meta$tf))
    }
    cat(paste0(
        'A Network object\n', 'with ', n_tfs, ' TFs and ',
        n_genes, ' target genes'
    ))
}

setMethod('show', 'Network', function(object) print(object))


#' Print Modules objects
#'
#' @rdname print
#' @export
#' @method print Modules
print.Modules <- function(object){
    n_mods <- length(object@features$genes_pos)
    cat(paste0(
        'An Modules object with ', n_mods, ' TF modules'
    ))
}

setMethod('show', 'Modules', function(object) print(object))


#' Print Regions objects
#'
#' @rdname print
#' @export
#' @method print Regions
print.Regions <- function(object){
    n_regs <- length(object@ranges)
    n_peaks <- length(unique(object@peaks))
    cat(paste0(
        'An Regions object\n', 'with ', n_regs, ' candidate genomic regions ',
        'in ', n_peaks, ' peaks'
    ))
}

setMethod('show', 'Regions', function(object) print(object))






