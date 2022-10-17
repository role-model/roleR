#' @title Parameters of one roleModel
#' @description An S4 class containing params for population sizes, rates of processes, the number of iterations
#' to run, and more 
#' 
#' @slot individuals_local number of individuals in local community
#' @slot individuals_meta number of individuals in meta community
#' @slot species_meta number of species in meta community
#' @slot speciation_local local speciation probability
#' @slot speciation_meta speciation rate in meta community
#' @slot extinction_meta extinction rate in meta community
#' @slot trait_sigma rate of Brownian trait evolution in the meta community
#' @slot env_sigma selectivity of environmental filter # might need another parameter for death rate in the face of neutrality
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
#' @details Params `init_type`, `niter`, `niterTimestep`, 
#'     `mutation_rate`,`equilib_escape`,and `num_basepairs` take a single value.
#'      All other params are numeric vectors containing either one value or `niter` values.
#'      If one value that value is used for all iterations of the model.
#'      If `niter`values a different sequential value is used for each iteration 
#'      
#' @examples 
#' Create a set of params
#' params <- roleParams(individuals_local = 100, individuals_meta = 1000,
#' species_meta = 10, speciation_local = 0.5, 
#' speciation_meta = 1, extinction_meta = 0.8, env_sigma = 0.5,
#' trait_sigma=1,comp_sigma = 0.1, dispersal_prob = 0.1, mutation_rate = 0.01,
#' equilib_escape = 1, num_basepairs = 250,
#' init_type = 'oceanic_island', niter = 2, niterTimestep = 2)
#' 
#' Use it to create a roleModel
#' model <- roleModel(params)
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
                           neut_delta = "numeric",
                           env_comp_delta = "numeric",
                           dispersal_prob = "numeric",
                           mutation_rate = "numeric" ,
                           equilib_escape = "numeric",
                           alpha = "numeric",
                           num_basepairs = "numeric",
                           init_type = "character", 
                           niter = 'integer', 
                           niterTimestep = 'integer'
                       )
                       )



# constructor
#' @rdname roleParams
#' @export

roleParams <- function(individuals_local,
                       individuals_meta,
                       species_meta,
                       speciation_local,
                       speciation_meta,
                       extinction_meta,
                       trait_sigma,
                       env_sigma,
                       comp_sigma,
                       neut_delta=NA,
                       env_comp_delta=NA,
                       dispersal_prob,
                       mutation_rate=NA,
                       equilib_escape=NA,
                       alpha=NA,
                       num_basepairs=NA,
                       init_type, 
                       niter, 
                       niterTimestep) {
    # set defaults - doesnt work right now 
    if(is.na(neut_delta)){neut_delta <- 0}
    if(is.na(env_comp_delta)){env_comp_delta <- 0.5}
    if(is.na(alpha)){alpha <- 1}
    if(is.na(equilib_escape)){equilib_escape <- 1}
    if(is.na(mutation_rate)){mutation_rate <- 0.001}
    if(is.na(num_basepairs)){num_basepairs <- 250}
    
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
            allParams[[i]] <- 0 #NA
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
            allParams$niterTimestep <- 10 * allParams$individuals_local[1]
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
               neut_delta = as.numeric(allParams$neut_delta),
               env_comp_delta = as.numeric(allParams$env_comp_delta),
               dispersal_prob = as.numeric(allParams$dispersal_prob),
               mutation_rate = as.numeric(allParams$mutation_rate),
               equilib_escape = as.numeric(allParams$equilib_escape),
               alpha = as.numeric(allParams$alpha),
               num_basepairs = as.numeric(allParams$num_basepairs),
               init_type = as.character(allParams$init_type), 
               niter = as.integer(allParams$niter),
               niterTimestep = as.integer(allParams$niterTimestep)))
}

# constructor
#' @rdname untbParams
#' @export

untbParams <- function(individuals_local,
                       individuals_meta,
                       species_meta,
                       speciation,
                       dispersal_prob,
                       init_type, 
                       niter, 
                       niterTimestep) {
    
    return(roleParams(
               individuals_local = individuals_local,
               individuals_meta = individuals_meta,
               species_meta = species_meta,
               speciation_local = speciation,
               speciation_meta = 0.8,
               extinction_meta = 0.05,
               trait_sigma = 1,
               env_sigma = 1,
               comp_sigma = 1,
               neut_delta = 1, # makes the model neutral by ignoring env and comp sigmas
               env_comp_delta = 0.5,
               dispersal_prob = dispersal_prob,
               mutation_rate = 0.01,
               equilib_escape = 1,
               alpha = 10,
               num_basepairs = 250,
               init_type = init_type, 
               niter = niter,
               niterTimestep = niterTimestep))
}