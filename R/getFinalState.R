#' @title Get the final state of a model run
#' @description gets the final state of a model from either a `roleExperiment` 
#'     or a `roleModel`
#' @param x the object to retrieve state from
#' @param modID optional argument to retrieve state for specific model(s) in a
#'     `roleExperiment` object; can be a vector of unique IDs
#' 
#' @details if `modID` is not supplied and `x` is a `roleExperiment` then a 
#'     list of `roleData` objects for each model run is returned
#' 
#' @export

setGeneric('getFinalState', 
           def = function(x, modID, ...) standardGeneric('getFinalState'), 
           signature = 'x')

setMethod('getFinalState', 
          signature = 'roleModel', 
          definition = function(x) {
              xlast <- x@modelSteps[length(x@modelSteps)]
              
              return(new('roleModel', 
                         params = x@params, 
                         modelSteps = xlast))
          }
)


setMethod('getFinalState', 
          signature = 'roleExperiment', 
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

