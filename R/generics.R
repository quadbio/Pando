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

#' @rdname GetGRN
#' @export GetGRN
GetGRN <- function(object, ...){
    UseMethod(generic = 'GetGRN', object = object)
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



