#' @title Get raw data from `roleData` objects.
#' @description Gets raw data according to the function specification.
#' @param x A `roleData` object.
#' @param ... Additional parameters, ignored
#' 
#' @return A list with one element. Inside that one element is a vector of the
#'     raw data.
#' 
#' @details These functions return a list of length one (with summary statistics as a vector
#'     inside that list) because the summary statistics themselves can be variable
#'     length. This is distinct from other summary statistic functions which will
#'     always return a fixed length value.  When using raw summary statistic
#'     functions in `getSumStats` they will be returned in a list column.
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
              list(tabulate(x@localComm@indSpecies))
          }
)

# abundance
setGeneric('rawSpAbundance', 
           def = function(x, ...) standardGeneric('rawSpAbundance'), 
           signature = 'x')

setMethod('rawSpAbundance', 
          signature = 'roleData', 
          definition = function(x) {
              list(tabulate(x@localComm@spAbund))
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
              list(x@localComm@indSpecies) #list(x@localComm@indSppTrt[, 1])
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
              list(x@localComm@indTrait)
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
              list(x@localComm@spGenDiv)  #list(x@localComm@sppGenDiv[, 1])
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
              list(x@localComm@indSeqs)
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

