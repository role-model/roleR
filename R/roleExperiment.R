#' @title roleExperiment - one or more models bundled collectively 
#' 
#' @description An S4 class to represent a self-enclosed modeling experiment 
#' 
#' 
#' @slot modelRuns a list of `roleData`objects
#' @slot allParams a list of `roleParams` to use for each model
# @slot context a named string vector that keep track of author metadata
#     It contains values for "author", "date", "description", "info", where 
#     each element is named by its respective string. When the model is saved 
#     with `writeRole` a text file is generated using this metadata
#' @slot info a data.frame with model parameter information for each saved 
#'     snapshot for all the models of the experiment 
#' @slot inits an optional list of initialized `roleModel` objects
#' 
#' 
#' @rdname roleExperiment-class
#' @include roleModel.R
#' @import methods
#' @export

setClass('roleExperiment',
         slots = c(info = 'data.frame',
                   modelRuns = 'list', 
                   allParams = 'list', 
                   inits = 'list'))


#' Create a RoLE experiment
#' @description Takes a list of parameter objects to initialize models for each
#'     parameter set.  The collection of models is seen as an experiment
#' 
#' @param allParams the list of model params to use for each model; can also be
#'     `rolePriors` object
#' @return a ready-to-run `roleExperiment`
#' 
#' @rdname roleExperiment
#' @import methods
#' @export

roleExperiment <- function(allParams) {
    # if given a `rolePriors` object, make it a list before proceeding 
    if(inherits(allParams, 'rolePriors')) {
        allParams <- list() # just a stub for now
    }
    
    # a list of `roleModel` objects
    allMods <- lapply(allParams, roleModel)
    
    # a list of `roleExperiment` objects converted from `allMods`
    allExps <- lapply(allMods, as, Class = 'roleExperiment')
    
    # combine list of `roleExperiments`
    outExp <- Reduce(rbind2, allExps)
    
    # write over the `inits` slot
    outExp@inits <- allMods
    
    return(outExp)
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
              origModID <- unique(x@info$mod_id)
              
              # update metadata and runs
              x@info <- x@info[i, j, ..., drop = FALSE]
              x@modelRuns <- x@modelRuns[i]
              
              # check to see if any unique model runs have been dropped
              theseDropped <- !(origModID %in% x@info$mod_id)
              
              # if so, drop them
              if(any(theseDropped)) {
                  i <- origModID[!theseDropped]
                  x@allParams <- x@allParams[i]
                  
                  # re-number IDs in metadata to remove any missing
                  x@info$mod_id <- order(x@info$mod_id)
              }
              
              return(x)
          })



#' Extract columns from `info` data.frame
#' @param x a `roleExperiment` object
#' @param name column name to extract
#' 
#' @name $
#' @aliases $,roleExperiment-method
#' @docType methods
#' @rdname roleExperiment-methods
#' 
setMethod("$", "roleExperiment", function(x, name) {
    if(name %in% names(x@info)) {
        return(x@info[, name])
    } else {
        stop('no $ method for object without metadata')
    }
})


# show method for `roleExperiment`
setMethod('show', signature = signature(object = 'roleExperiment'),
          definition = function(object) {
              nmod <- length(unique(object@info$mod_id))
              
              is_run <- !is.null(object@modelRuns[[2]]) # see if model has been run
              
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
              
              # keep track of growing mod_id max index
              j <- max(out@info$mod_id)
              
              # augment mod_id
              thisEx@info$mod_id <- thisEx@info$mod_id + j
              
              # combine with `out`
              out@info <- rbind(out@info, thisEx@info)
              out@modelRuns <- c(out@modelRuns, thisEx@modelRuns)
              out@allParams <- c(out@allParams, thisEx@allParams)
              out@inits <- c(out@inits, thisEx@inits)
              
              return(out)
          }
)




setMethod('rbind2', signature = c('roleExperiment', 'missing'), 
          definition = function(x, y) {
              return(x)
          })


