#' @title Get summary statistics for RoLE objects.
#' @description Applies different summary statistics functions to `roleExperiment`,
#'     `roleModel`, or `roleData` objects.
#' @param x The object to calculate summary statistics on. 
#' @param funs A named list of functions to calculate the summary statistics; can be a 
#'     named list with a single function or many functions, but must be a 
#'     named list of functions. If unspecified, defaults to including several standard stats. 
#'     Included diversity index functions in roleR are: `hillAbund`, 
#'     `hillGenetic`,`hillPhylo` & `richness`. Included raw stats functions are `rawAbundance`,
#'     `rawSpAbundance`,`rawSppId`,`rawTraits`,`rawGenDiv`,`rawSeqs`,
#'     `rawBranchLengths` & `rawApePhylo`.
#'     
#' @param moreArgs An optional named list of additional arguments to pass to the
#'     functions listed in `funs`. If given, list names must match those in 
#'     `funs`. Note: not all names in `funs` need to appear in `moreArgs`
#' @param ... ignored
#' 
#' @details Users can define their own functions, so long as they work on any
#'     object of class `roleData`. 
#' 
#' @return data.frame containing summary stats where each row is a model snapshot
#' and each column is a summary stat requested by a function provided to `funs`
#' 
#' @include roleData.R
#' 
#' @examples 
#' # get the species richness 
#' # rich_stat <- getSumStats(model, funs = list(rich = richness))
#' # get many default summary stats
#' # stats <- getSumStats(model)
#' 
#' @rdname getSumStats
#' @export

setGeneric('getSumStats', 
           def = function(x, funs=.getDefaultSumStatsFuns(), moreArgs, ...) standardGeneric('getSumStats'), 
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
              o <- do.call(rbind, o)
             
              o$iteration <- seq(from=0,by=x@params@niterTimestep,length.out = (x@params@niter / x@params@niterTimestep) + 1)
              return(o)
          }
)


# method for roleExperiment
#' @rdname getSumStats
#' @export

setMethod('getSumStats', 
          signature = 'roleExperiment', 
          definition = function(x, funs, moreArgs) {
              
              ss <- data.frame()
              for(r in 1:length(x@modelRuns)){
                  s <- getSumStats(x@modelRuns[[r]],funs)
                  s$run_num <- r
                  ss <- rbind(ss,s)
              }
              return(ss)
          }
)

#' Get SumStats mean
#'
#' @param x x
#' @param funs funs
#'
#' @return something
getSumStatsMean <- function(x, funs){
    ss <- list()
    for(r in 1:length(x@modelRuns)){
        s <- getSumStats(x@modelRuns[[r]],funs)
        s$run_num <- r
        ss <- append(ss,list(s))
    }
}

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

# helper function that returns the default list of summary stats and names
# used in the constructor of getSumStats
.getDefaultSumStatsFuns <- function(){
    return(list(hill_abund=hillAbund, hill_gen=hillGenetic, hill_trait = hillTrait, hill_phy = hillPhylo,richness=richness
                ))
    
}