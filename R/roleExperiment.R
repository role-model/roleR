#' @title roleExperiment - summaries of one or more models bundled collectively 
#' 
#' @description An S4 class to represent a self-enclosed modeling experiment 
#' 
#' 
#' @slot allParams a list of `roleParams` objects used for each model
#' @slot results a data.frame with one row per saved snapshot for all the 
#'     models of the experiment. The first columns are `model_id` and the model
#'     information from the `info` slot of each initial `roleModel`; once the
#'     experiment has been run, the summary statistics are added as columns
#' @slot funs a named list of summary stat functions
#' @slot inits a list of initialized `roleModel` objects used as
#'     initial conditions (by default, calculated internally, but can be 
#'     passed as an optional argument)
#' 
#' 
#' @rdname roleExperiment-class
#' @include roleModel.R
#' @import methods
#' @export

setClass('roleExperiment',
         slots = c(results = 'data.frame',
                   allParams = 'list', 
                   funs = "list",
                   inits = 'list'))


#' Create a RoLE experiment
#' @description Takes a list of parameter objects to build models for each
#'     parameter set. The collection of models is seen as an experiment. Raw
#'     `roleModel` objects will not be returned, rather summaries (indicated by
#'     `funs`) will be returned
#' 
#' @param allParams the list of model params to use for each model; can also be
#'     `rolePriors` object
#' @param funs named list of summary stat functions
#' @param inits optional list of initialized `roleModel` objects used as 
#'     initial conditions
#' @return a ready-to-run `roleExperiment`
#' 
#' @rdname roleExperiment
#' @import methods
#' @export

roleExperiment <- function(allParams, funs = NULL, 
                           inits = NULL) {
    # if given a `rolePriors` object, make it a list before proceeding 
    if(inherits(allParams, 'rolePriors')) {
        allParams <- list() # just a stub for now
    }
    
    if(length(allParams) == 0) stop("`allParams` must not be empty")
    
    if(is.null(funs)) funs <- .getDefaultSumStatsFuns()
    
    # a list of `roleModel` objects
    if(is.null(inits)) {
        inits <- lapply(allParams, roleModel)
    } else {
        # make sure params are correct for user-provided `roleModel` objects
        for(i in seq_along(allParams)) {
            inits[[i]]@params <- allParams[[i]]
        }
    }
    
    # extract `info` slot from all `roleModel` objects
    allResults <- lapply(seq_along(inits), function(i) {
        d <- cbind(model_id = i, slot(inits[[i]], name = "info"))
        d
    }) |> 
        do.call(rbind, args = _)
    
    return(new("roleExperiment", 
               allParams =  allParams, 
               results = allResults,
               funs = funs,
               inits = inits))
}




#' Extract parts of `roleExperiment` object
#' @param x a `roleExperiment` object
#' @param i row index
#' @param j column index
#' @param drop boolean, whether to drop dimensions
#' @param ... additional arguemnts 
#' 
# @name [
#' @aliases [,roleExperiment-method
#' @docType methods
#' @rdname roleExperiment-methods
#'
setMethod("[", 
          signature(x = "roleExperiment", i = "ANY", j = "ANY"),
          function(x, i, j, ..., drop = FALSE) {
              iMissing <- missing(i)
              jMissing <- missing(j)
              nargs <- nargs() # e.g., a[3,] gives 2 for nargs, a[3] gives 1.
              
              if(iMissing && jMissing) {
                  i <- TRUE
                  j <- TRUE
              } else if(jMissing && !iMissing) { 
                  if (nargs == 2) {
                      j <- i
                      i <- TRUE
                  } else {
                      j <- TRUE
                  }
              } else if(iMissing && !jMissing)
                  i <- TRUE
              
              if(is.matrix(i)) {
                  stop('matrix argument not supported')
              }
              
              if(any(is.na(i))) {
                  stop('NAs not permitted in row index')
              }
              
              # vector of unique model run IDs before selecting
              origModID <- unique(x@results$model_id)
              
              # update results
              x@results <- x@results[i, j, ..., drop = FALSE]
              
              # check to see if any unique models have been dropped
              theseDropped <- !(origModID %in% x@results$model_id)
              
              # if so, drop them
              if(any(theseDropped)) {
                  i <- origModID[!theseDropped]
                  x@allParams <- x@allParams[i]
                  x@inits <- x@inits[i]
                  
                  # re-number IDs in `results` to remove any missing
                  x@results$model_id <- order(x@results$model_id)
              }
              
              return(x)
          })



#' Extract columns from `results` data.frame
#' @param x a `roleExperiment` object
#' @param name column name to extract
#' 
#' @name $
#' @aliases $,roleExperiment-method
#' @docType methods
#' @rdname roleExperiment-methods
#' 
setMethod("$", "roleExperiment", function(x, name) {
    if(name %in% names(x@results)) {
        return(x@results[, name])
    } else {
        stop('no column named ', name, ' in results')
    }
})


# show method for `roleExperiment`
setMethod('show', signature = signature(object = 'roleExperiment'),
          definition = function(object) {
              nmod <- length(unique(object@results$model_id))
              
              # the experiment has been run if `results` has columns beyond
              # `model_id` and the model information
              is_run <- ncol(object@results) > 
                  1 + ncol(object@inits[[1]]@info)
              
              if(is_run){
                  run_str <- "completed (run)"
              }
              
              else{
                  run_str <- "not-yet-run"
              }
              
              cat(sprintf('%s RoLE experiment with %s unique model%s',
                          run_str, 
                          nmod, 
                          ifelse(nmod == 1, '', 's')), 
                  '\n')
          }
)



# rbind method for `roleExperiment` class

setMethod('rbind2', signature = c('roleExperiment', 'roleExperiment'), 
          definition = function(x, y) {
              out <- x
              thisEx <- y
              
              # keep track of growing model_id max index
              j <- max(out@results$model_id)
              
              # augment model_id
              thisEx@results$model_id <- thisEx@results$model_id + j
              
              # combine with `out`
              out@results <- rbind(out@results, thisEx@results)
              out@allParams <- c(out@allParams, thisEx@allParams)
              out@inits <- c(out@inits, thisEx@inits)
              
              return(out)
          }
)




setMethod('rbind2', signature = c('roleExperiment', 'missing'), 
          definition = function(x, y) {
              return(x)
          })


