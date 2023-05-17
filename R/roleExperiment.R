#' @title roleExperiment - one or more models bundled collectively 
#' 
#' @description An S4 class to represent a self-enclosed modeling experiment 
#' 
<<<<<<< HEAD
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
=======
#' @slot allFuns a list, something
#' @slot modelRuns a list of `roleModel`objects
#' @slot allParams a list of `roleParams` to use for each model
#' @slot context a named string vector that keep track of author metadata
#' It contains values for "author", "date", "description", "info", where each element is named by its respective string.
#' When the model is saved with `writeRole` a text file is generated using this metadata
#' @slot info data.frame of summarizing metadata for all the models of the experiment 
#' (Need to chat with Andy about the exact intentions of this before writing more here)
#' @include roleModel.R roleParams.R
>>>>>>> main
#' @examples 
#' # create and run a roleExperiment with three different levels of dispersal 
# p1 <- roleParams(dispersal_prob = 0.1)
# p2 <- roleParams(dispersal_prob = 0.2)
# p3 <- roleParams(dispersal_prob = 0.3)
# expr <- roleExperiment(list(p1, p2, p3))
# expr <- runRole(expr)
#' 
#' @rdname roleExperiment
#' @export

setClass('roleExperiment',
         slots = c(info = 'data.frame',
                   modelRuns = 'list', 
                   allParams = 'list', 
                   inits = 'list'))


#' @param allParams the list of model params to use for each model; can also be
#'     `rolePriors` object
#' @return a ready-to-run `roleExperiment`
#' 
#' @rdname roleExperiment
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


# set coercion method from `roleModel` to `roleExperiment`
setAs(from = 'roleModel', to = 'roleExperiment',
      def = function(from) {
          # get metadata
          pout <- from@info
          
          # put an index column in that data.frame for model 
          pout <- cbind(mod_id = 1, pout)
          
          return(new('roleExperiment', 
                     info = pout,
                     modelRuns = from@modelSteps, 
                     allParams = list(from@params), 
                     inits = list()))
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


<<<<<<< HEAD
=======
#' Repeat an S4 object
#' @description temporary helper to rep an S4 object
#' used to create roleExperiments with many models of the same params
#' i.e. roleExperiment(repS4(p,100)) does this 100 times using roleParams p
#' @param obj an s4
#' @param n times
#'
#' @return obj repeated n times
#' @export
#'
repS4 <- function(obj, n) {
    # Create an empty list to store replicated objects
    objs <- vector("list", n)
    
    # Replicate the object and store in the list
    for (i in 1:n) {
        objs[[i]] <- obj
    }
    
    # Return the list of replicated objects
    return(objs)
}

#' @include rolePhylo.R
# set coercion method from ape::phylo to roleR::rolePhylo
setAs(from = 'phylo', to = 'rolePhylo',
      def = function(from) {
          # extract number of times
          n <- ape::Ntip(from)
          
          # extract edge matrix and edge lengths
          e <- from$edge
          l <- from$edge.length
          
          # extract tip labels
          tipNames <- from$tip.label
          
          # calculate alive or not
          tipAge <- ape::node.depth.edgelength(from)[1:n]
          
          alive <- rep(TRUE, n)
          alive[tipAge < max(tipAge)] <- FALSE
          
          # set default scale
          scale <- 1
          
          return(rolePhylo(n = n, e = e, l = l, alive = alive,
                           tipNames = tipNames, scale = scale))
      }
)


# set coercion method from roleR::rolePhylo to ape::phylo
setAs(from = 'rolePhylo', to = 'phylo',
      def = function(from) {
          i <- 2 * (from@n - 1)
          
          y <- list(edge = from@e[1:i, ], edge.length = from@l[1:i],
                    tip.label = from@tipNames[1:from@n],
                    Nnode = from@n - 1)
          
          # make any possible 0 or negative edge lengths equal to
          # very small number
          y$edge.length[y$edge.length <= 0] <- .Machine$double.eps
          
          class(y) <- 'phylo'
          
          return(y)
      }
)
>>>>>>> main
