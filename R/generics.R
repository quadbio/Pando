#' @rdname initiate_grn
#' @export initiate_grn
initiate_grn <- function(object, ...){
    UseMethod(generic = 'initiate_grn', object = object)
}

#' @rdname find_motifs
#' @export find_motifs
find_motifs <- function(object, ...){
    UseMethod(generic = 'find_motifs', object = object)
}

#' @rdname infer_grn
#' @export infer_grn
infer_grn <- function(object, ...){
    UseMethod(generic = 'infer_grn', object = object)
}

#' @rdname fit_grn_models
#' @export fit_grn_models
fit_grn_models <- function(object, ...){
    UseMethod(generic = 'fit_grn_models', object = object)
}

#' @rdname plot_gof
#' @export plot_gof
plot_gof <- function(object, ...){
    UseMethod(generic = 'plot_gof', object = object)
}

#' @rdname plot_module_metrics
#' @export plot_module_metrics
plot_module_metrics <- function(object, ...){
    UseMethod(generic = 'plot_module_metrics', object = object)
}

#' @rdname get_network_graph
#' @export get_network_graph
get_network_graph <- function(object, ...){
    UseMethod(generic = 'get_network_graph', object = object)
}

#' @rdname plot_network_graph
#' @export plot_network_graph
plot_network_graph <- function(object, ...){
    UseMethod(generic = 'plot_network_graph', object = object)
}

#' @rdname get_tf_network
#' @export get_tf_network
get_tf_network <- function(object, ...){
    UseMethod(generic = 'get_tf_network', object = object)
}

#' @rdname plot_tf_network
#' @export plot_tf_network
plot_tf_network <- function(object, ...){
    UseMethod(generic = 'plot_tf_network', object = object)
}

#' @rdname NetworkGraph
#' @export NetworkGraph
NetworkGraph <- function(object, ...){
    UseMethod(generic = 'NetworkGraph', object = object)
}

#' @rdname GetGRN
#' @export GetGRN
GetGRN <- function(object, ...){
    UseMethod(generic = 'GetGRN', object = object)
}

#' @rdname GetNetwork
#' @export GetNetwork
GetNetwork <- function(object, ...){
    UseMethod(generic = 'GetNetwork', object = object)
}

#' @rdname NetworkFeatures
#' @export NetworkFeatures
NetworkFeatures <- function(object, ...){
    UseMethod(generic = 'NetworkFeatures', object = object)
}

#' @rdname NetworkRegions
#' @export NetworkRegions
NetworkRegions <- function(object, ...){
    UseMethod(generic = 'NetworkRegions', object = object)
}

#' @rdname Params
#' @export Params
Params <- function(object, ...){
    UseMethod(generic = 'Params', object = object)
}

#' @rdname NetworkParams
#' @export NetworkParams
NetworkParams <- function(object, ...){
    UseMethod(generic = 'NetworkParams', object = object)
}

#' @rdname NetworkTFs
#' @export NetworkTFs
NetworkTFs <- function(object, ...){
    UseMethod(generic = 'NetworkTFs', object = object)
}

#' @rdname NetworkModules
#' @export NetworkModules
NetworkModules <- function(object, ...){
    UseMethod(generic = 'NetworkModules', object = object)
}

#' @rdname DefaultNetwork
#' @export DefaultNetwork
DefaultNetwork <- function(object, ...){
    UseMethod(generic = 'DefaultNetwork', object = object)
}

#' @rdname gof
#' @export gof
gof <- function(object, ...){
    UseMethod(generic = 'gof', object = object)
}

#' @rdname find_modules
#' @export find_modules
find_modules <- function(object, ...){
    UseMethod(generic = 'find_modules', object = object)
}

#' @rdname GetAssaySummary
#' @export GetAssaySummary
GetAssaySummary <- function(object, ...){
    UseMethod(generic = 'GetAssaySummary', object = object)
}



