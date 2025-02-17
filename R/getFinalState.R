#' @title Get the final state of a model run.
#' @description Gets the final state of a model from either a `roleExperiment` 
#'     or a `roleModel`.
#' @param x The `roleExperiment` or `roleModel` object to retrieve the final state from.
#' @param modID Optional argument to retrieve the final state for specific model(s) in a
#'     `roleExperiment` object; can be a vector of unique IDs. 
#' @param ... passed
#' 
#' @details If `modID` is not supplied and `x` is a `roleExperiment`, then a 
#'     list of `roleData` objects for each model run is returned.
#' 
#' @include roleModel.R roleExperiment.R
#' @rdname getFinalState
#' @export

setGeneric('getFinalState', 
           def = function(x, modID, ...) standardGeneric('getFinalState'), 
           signature = 'x')


#' Get final state of roleModel
#' @name getFinalState
#' @aliases getFinalState,roleModel-method
#' @docType methods
#' @rdname getFinalState
#' @export

setMethod('getFinalState', 
          signature(x = 'roleModel'), 
          definition = function(x) {
              xlast <- x@modelSteps[length(x@modelSteps)]
              
              return(new('roleModel', 
                         params = x@params, 
                         modelSteps = xlast))
          }
)


#' Get final state of roleExperiment
#' @name getFinalState
#' @aliases getFinalState,roleExperiment-method
#' @docType methods
#' @rdname getFinalState
setMethod('getFinalState', 
          signature(x ='roleExperiment'), 
          definition = function(x, modID) {
              # if IDs are given, subset based on those
              if(!missing(modID)) {
                  x <- x[x@experimentMeta$mod_id %in% modID, ]
              } else {
                  modID <- unique(x@experimentMeta$mod_id)
              }
              
              # calculate index of final time step for each model
              iFinal <- sapply(modID, function(i) {
                  j <- x@experimentMeta$iterations[x@experimentMeta$mod_id == i]
                  jMax <- max(j)
                  
                  k <- x@experimentMeta$iterations == jMax & 
                      x@experimentMeta$mod_id == i
                  
                  return(which(k))
              })
              
              return(x[iFinal, ])
          }
)

