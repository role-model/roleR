#' @title Get raw data from `roleData` objects
#' @description Gets raw data
#' @param x the `roleData` object to get raw data from
#' @param ... additional parameters, ignored
#' 
#' @return a list with one element, inside that one element is a vector of the
#'     raw data
#' 
#' @details These functions return a list of length one (with summary statistics 
#'     inside that list) because the summary stats themselves could be variable
#'     length. This is distinct from other summary stat functions which will
#'     always return a fixed length value.  When using raw summary stat 
#'     functions in `getSumStats` they will be returned in a list column
#' 
#' @rdname raw-sumStats
#' @export

# abundance
setGeneric('rawAbundance', 
           def = function(x, ...) standardGeneric('rawAbundance'), 
           signature = 'x')

setMethod('rawAbundance', 
          signature = 'roleData', 
          definition = function(x) {
              list(tabulate(x@localComm@indSppTrt[, 1]))
          }
)


# spp IDs
#' @rdname raw-sumStats
#' @export

setGeneric('rawSppID', 
           def = function(x, ...) standardGeneric('rawSppID'), 
           signature = 'x')

setMethod('rawSppID', 
          signature = 'roleData', 
          definition = function(x) {
              list(x@localComm@indSppTrt[, 1])
          }
)

# traits
#' @rdname raw-sumStats
#' @export

setGeneric('rawTraits', 
           def = function(x, ...) standardGeneric('rawTraits'), 
           signature = 'x')

setMethod('rawTraits', 
          signature = 'roleData', 
          definition = function(x) {
              list(x@localComm@indSppTrt[, 2])
          }
)


# gen div
#' @rdname raw-sumStats
#' @export

setGeneric('rawGenDiv', 
           def = function(x, ...) standardGeneric('rawGenDiv'), 
           signature = 'x')

setMethod('rawGenDiv', 
          signature = 'roleData', 
          definition = function(x) {
              list(x@localComm@sppGenDiv[, 1])
          }
)



# gen seqs
#' @rdname raw-sumStats
#' @export

setGeneric('rawSeqs', 
           def = function(x, ...) standardGeneric('rawSeqs'), 
           signature = 'x')

setMethod('rawSeqs', 
          signature = 'roleData', 
          definition = function(x) {
              list(x@localComm@indSeqs[, 1])
          }
)


# branch lengths
#' @rdname raw-sumStats
#' @export

setGeneric('rawBranchLengths', 
           def = function(x, ...) standardGeneric('rawBranchLengths'), 
           signature = 'x')

setMethod('rawBranchLengths', 
          signature = 'roleData', 
          definition = function(x) {
              list(x@phylo@l[x@phylo@e[, 1] != -1])
          }
)


# ape phylo
#' @rdname raw-sumStats
#' @export

setGeneric('rawApePhylo', 
           def = function(x, ...) standardGeneric('rawApePhylo'), 
           signature = 'x')

setMethod('rawApePhylo', 
          signature = 'roleData', 
          definition = function(x) {
              list(as(x@phylo, 'phylo'))
          }
)
