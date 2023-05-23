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
#' @examples 
#' # create and run a roleExperiment with three different levels of dispersal 
# p1 <- roleParams(dispersal_prob = 0.1)
# p2 <- roleParams(dispersal_prob = 0.2)
# p3 <- roleParams(dispersal_prob = 0.3)
# expr <- roleExperiment(list(p1, p2, p3))
# expr <- runRole(expr)
#' 
#' @rdname roleExperiment
#' @include roleModel.R
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
