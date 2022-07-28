#' @title Priors of one roleModel
#' 
#' @slot individuals_local number of individuals in local community
#' @slot individuals_meta number of individuals in meta community
#' @slot species_meta number of species in meta community
#' @slot speciation_local local speciation probability
#' @slot speciation_meta speciation rate in meta community
#' @slot extinction_meta extinction rate in meta community
#' @slot trait_sigma rate of Brownian trait evolution in the meta community
#' @slot env_sigma selectivity of environmental filter
#' @slot comp_sigma selectivity of competition
#' @slot dispersal_prob local dispersal probability
#' 
#' @slot mutation_rate mutation rate
#' @slot equilib_escape proportion of equilibrium achieved
#' @slot num_basepairs number of basepairs
#' 
#' @slot init_type initialization routine; a single character string either 
#'     "oceanic_island" or "bridge_island"
#'     
#' @slot niter an integer specifying the number of iterations 
#' @slot niterTimestep an integer specifying the frequency (in numbers of 
#'     iterations) at which the model state is saved
#' 
#' @details Each parameter other than `init_type`, `niter`, and `niterTimestep`
#'     is a vector containing either one value or `niter` values; `init_type`, 
#'     `niter`, and `niterTimestep` can only take a single value (others can
#'     only take single vals too...clear this up)
#' 
#' @rdname roleIterFuns
#' @export

roleIterFuns <- setClass('roleIterFuns',
                       slots = c(
                           individuals_local = "list",
                           individuals_meta = "list",
                           species_meta = "list",
                           speciation_local = "list",
                           speciation_meta = "list",
                           extinction_meta = "list",
                           trait_sigma = "list",
                           env_sigma = "list",
                           comp_sigma = "list",
                           dispersal_prob = "list",
                           mutation_rate = "list",
                           equilib_escape = "list",
                           num_basepairs = "list",
                           init_type = "list", 
                           niter = 'list', 
                           niterTimestep = 'list'
                       ))



# constructor
#' @rdname roleIterFuns
#' @export

roleIterFuns <- function(individuals_local,
                       individuals_meta,
                       species_meta,
                       speciation_local,
                       speciation_meta,
                       extinction_meta,
                       trait_sigma,
                       env_sigma,
                       comp_sigma,
                       dispersal_prob,
                       mutation_rate,
                       equilib_escape,
                       num_basepairs,
                       init_type, 
                       niter, 
                       niterTimestep) {
    # check that `niter` is given correctly
    if(missing(niter) | length(niter) > 1) {
        stop('must supply a single value for `niter`')
    }
    
    # check for missing params
    allParams <- as.list(environment())
    withVal <- names(as.list(match.call())[-1])
    noVal <- names(allParams[!(names(allParams) %in% withVal)])
    
    # loop over params:
    # those that need to be `niter` long, make them that
    # those that are missing, make them `NA`
    singleValParams <- c('individuals_meta', 'species_meta',
                         'speciation_meta', 'extinction_meta', 'trait_sigma',
                         'equilib_escape', 'num_basepairs', 'init_type', 
                         'niter', 'niterTimestep')
    
    for(i in 1:length(allParams)) {
        if(names(allParams[i]) %in% noVal) {
            allParams[[i]] <- NA
        }
        
        if(names(allParams[i]) %in% singleValParams) {
            if(length(allParams[[i]]) > 1) {
                stop(sprintf('`%s` must be a single value'))
            } 
        } else {
            if(length(allParams[[i]]) == 1) {
                allParams[[i]] <- rep(allParams[[i]], niter)
            } else if(length(allParams[[i]]) < niter) {
                stop(sprintf('`%s` must be a single value', names(allParams[i])), 
                     'or a vector of values exactly `niter` long')
            }
        }
    }
    
    # if `niterTimestep` was missing, make it a meaningful multiple of `niter`
    if(is.na(allParams$niterTimestep)) {
        if(allParams$niter < 10 * allParams$individuals_local[1]) {
            allParams$niterTimestep <- allParams$niter
        } else {
            allParams$niterTimestep <- 10 * allParams$individuals_local
        }
    }
    
    return(new('roleIterFuns', 
               individuals_local = allParams$individuals_local,
               individuals_meta = allParams$individuals_meta,
               species_meta = allParams$species_meta,
               speciation_local = (allParams$speciation_local),
               speciation_meta = (allParams$speciation_meta),
               extinction_meta = (allParams$extinction_meta),
               trait_sigma = (allParams$trait_sigma),
               env_sigma = (allParams$env_sigma),
               comp_sigma = (allParams$comp_sigma),
               dispersal_prob = (allParams$dispersal_prob),
               mutation_rate = (allParams$mutation_rate),
               equilib_escape = (allParams$equilib_escape),
               num_basepairs = (allParams$num_basepairs),
               init_type = allParams$init_type, 
               niter = allParams$niter,
               niterTimestep = allParams$niterTimestep))
}

# todo add checks that functions take an int iter