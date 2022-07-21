#' @title Get summary statistics for RoLE objects
#' @description Applies different summary stats functions to `roleExperiment`,
#'     `roleModel`, or `roleData`
#' @param x the object to calculate sum stats across
#' @param funs a named list of function to calculate the sum stats; can be a 
#'     named list with a single function or many functions, but must be a 
#'     named list of functions
#' @param moreArgs an optional named list of additional arguments to pass to the
#'     functions listed in `funs`; if given, list names must match those in 
#'     `funs`. Note: not all names in `funs` need to appear in `moreArgs`
#' @param ... ignored
#' 
#' @details users can define their own functions, so long as they work on any
#'     object of class `roleData`
#' 
#' @rdname getSumStats
#' @export

setGeneric('getSumStats', 
           def = function(x, funs, moreArgs, ...) standardGeneric('getSumStats'), 
           signature = 'x')


# method for roleData
#' @rdname getSumStats
#' @export

setMethod('getSumStats', 
          signature = 'roleData', 
          definition = function(x, funs, moreArgs) {
              bigFun <- .funApply(funs, moreArgs)
              
              return(bigFun(x))
          }
)


# method for roleModel
#' @rdname getSumStats
#' @export

setMethod('getSumStats', 
          signature = 'roleModel', 
          definition = function(x, funs, moreArgs) {
              bigFun <- .funApply(funs, moreArgs)
              
              o <- lapply(x@modelSteps, bigFun)
              
              return(do.call(rbind, o))
          }
)


# method for roleExperiment
#' @rdname getSumStats
#' @export

setMethod('getSumStats', 
          signature = 'roleExperiment', 
          definition = function(x, funs, moreArgs) {
              bigFun <- .funApply(funs, moreArgs)
              
              o <- lapply(x@modelRuns, bigFun)
              
              return(do.call(rbind, o))
          }
)



# helper function to make a new synthetic function from a list of functions 
.funApply <- function(funs, moreArgs) {
    if(missing(moreArgs)) moreArgs <- list(NULL)
    
    newFun <- function(x) {
        o <- lapply(1:length(funs), function(i) {
            # get function from list
            f <- funs[[i]]
            fname <- names(funs[i])
            
            # get more args (if they exist)
            fargs <- moreArgs[[fname]]
            if(is.null(fargs)) {
                allArgs <- list(x = x)
            } else {
                allArgs <- c(list(x = x), fargs)
            }
            
            # browser()
            # function value
            val <- do.call(f, allArgs)
            
            # name the values
            if(length(val) > 1) {
                ename <- names(val)
                
                if(is.null(ename)) {
                    ename <- 1:length(val)
                }
                
                ename <- paste(fname, 1:length(val), sep = '_')
            } else {
                ename <- fname
            }
            
            names(val) <- ename
            
            # return as a data.frame
            if('list' %in% class(val)) {
                # if val is a list, must protect it to make a list column
                df <- data.frame(I(val))
                names(df) <- names(val)
                return(df)
            } else {
                return(as.data.frame(as.list(val)))
            }
        })
        
        # combine data.frames and return
        o <- do.call(cbind, o)
        rownames(o) <- NULL
        
        return(o)
    }
    
    # the actual final return of `.funApply` is a function that runs the 
    # above code
    return(newFun)
}

