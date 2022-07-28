#' @title Parameters of one roleModel
#' @description An S4 class for roleModel parameters
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
#' @rdname roleParams
#' @export

roleParams <- setClass('roleParams',
                       slots = c(
                           individuals_local = "numeric",
                           individuals_meta = "numeric",
                           species_meta = "numeric",
                           speciation_local = "numeric",
                           speciation_meta = "numeric",
                           extinction_meta = "numeric",
                           trait_sigma = "numeric",
                           env_sigma = "numeric",
                           comp_sigma = "numeric",
                           dispersal_prob = "numeric",
                           mutation_rate = "numeric" ,
                           equilib_escape = "numeric",
                           alpha = "numeric",
                           num_basepairs = "numeric",
                           init_type = "character", 
                           niter = 'integer', 
                           niterTimestep = 'integer'
                       ))



# constructor
#' @rdname roleParams
#' @export

roleParams <- function(individuals_local=100,
                       individuals_meta=1000,
                       species_meta=500,
                       speciation_local=0.1,
                       speciation_meta=0.5,
                       extinction_meta=0.5,
                       trait_sigma=1,
                       env_sigma=0.1,
                       comp_sigma=0.1,
                       dispersal_prob=0.1,
                       mutation_rate=0.01,
                       equilib_escape=1,
                       num_basepairs=250,
                       init_type='oceanic_island', 
                       niter=100, 
                       niterTimestep=10) {
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
    
    return(new('roleParams', 
               individuals_local = as.numeric(allParams$individuals_local),
               individuals_meta = as.numeric(allParams$individuals_meta),
               species_meta = as.numeric(allParams$species_meta),
               speciation_local = as.numeric(allParams$speciation_local),
               speciation_meta = as.numeric(allParams$speciation_meta),
               extinction_meta = as.numeric(allParams$extinction_meta),
               trait_sigma = as.numeric(allParams$trait_sigma),
               env_sigma = as.numeric(allParams$env_sigma),
               comp_sigma = as.numeric(allParams$comp_sigma),
               dispersal_prob = as.numeric(allParams$dispersal_prob),
               mutation_rate = as.numeric(allParams$mutation_rate),
               equilib_escape = as.numeric(allParams$equilib_escape),
               num_basepairs = as.numeric(allParams$num_basepairs),
               init_type = as.character(allParams$init_type), 
               niter = as.integer(allParams$niter),
               niterTimestep = as.integer(allParams$niterTimestep)))
}
