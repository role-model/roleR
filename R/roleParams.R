#' @title Parameters of one roleModel
#' @description An S4 class containing params for population sizes, rates of processes, the number of iterations
#' to run, and much more. 
#' 
#' @slot individuals_local Number of individuals in the local community (J).
#' Determines the length of raw individual-level trait, abundance, and genetic data vectors returned for simulations.
#' When the model is initialized, all individuals are of a single species from the metacommunity.
#' @slot individuals_meta Number of individuals in the metacommunity.
#' Used in generating a log series of the initial abundances of species in the metacommunity. 
#' @slot species_meta Number of species in the metacommunity. 
#' Determines the initial size of the phylogeny.
#' 
#' @slot speciation_local Probability of speciation occurring in the local community at each time step.
#' The new local species can be either from a birth in the local community or an immigration from the metacommunity -
#' the `dispersal prob` param determines the chance of each of these scenarios.
#' @slot speciation_meta rate of speciation in meta community.
#' This is used during the simulation of the start phylogeny before the model is run
#' The sum of speciation_meta and extinction_meta is the average lifetime of phylo branches, 
#' and the larger this value the less new individual traits will deviate from
#' the parent's individual trait (if parent is in the local comm) or the metacommunity mean (if parent in the metacomm)
#' @slot extinction_meta Rate of extinction in the metacommunity.
#' Like speciation_meta, used in the initial phylogeny simulation and in relation to trait deviation.
#' @slot dispersal_prob Probability of dispersal (immigration) occurring from the metacommunity to the local local community at each time step.
#' Every time step, either birth or immigration happens, so the probability of birth is 1 minus the dispersal_prob.
#' @slot env_optim 
#' @slot trait_sigma Rate of Brownian trait evolution in the metacommunity. 
#' Determines how much the trait of a new individual deviates from its parent, i.e. how fast traits change.
#' @slot env_sigma Selectivity of the environmental filter, i.e. how strongly the environment selects which trait values are a good match for it.
#' The larger the value, the less chance an individual will survive if it's far from the trait optimum (which is 0).
#' @slot comp_sigma the size of the competition kernel: how strongly distance from the traits of other species 
#' selects which trait values are likely to survive.
#' The larger the value, the greater the chance an individual will survive if it is close in trait space to others
#' @slot neut_delta the degree of neutrality - if 1 there is no environmental filtering or competition, 
#' if 0 there is full competition and filtering.
#' @slot env_comp_delta a slider between environmental filtering and competition - if 1 there is full filtering and no competition,
#' if 0 there is full competition and no filtering.
#' 
#' @slot mutation_rate Rate of sequence mutation to use in genetic simulations.
#' @slot equilib_escape Proportion of equilibrium required to halt the model as it is running and return it
#' @slot alpha alpha parameter
#' @slot num_basepairs Number of basepairs to use in genetic simulations. Genetic simulations are currently single-locus.
#' 
#' @slot init_type The biological model used to initialize; a single character string that can be either "oceanic_island", "bridge_island", or "bare_island."
#' The bridge island model has the initial individuals in the local community arriving through a land bridge, while the oceanic has no bridge and is populated by a single dispersal event.
#' Thus, for oceanic island models, all individuals are of a SINGLE species sampled proportional to meta community species abundances, 
#' while in bridge island models, individuals are sampled of MANY species proportional to their abundances.
#' The bare island model is related to the oceanic island model, but instead of starting with "individuals_local" individuals of the sole sampled species, only 1 individual 
#' of that species appears and the rest of the space is filled with placeholder "rocks" representing unfilled space.
#' @slot niter An integer specifying the number of time steps for the model to run.
#' @slot niterTimestep An integer specifying the frequency (in numbers of 
#'     iterations) at which the model state is snapshot and saved in a model's model steps object.
#' 
#' @details Params `init_type`, `niter`, `niterTimestep`, 
#'     `mutation_rate`,`equilib_escape`,and `num_basepairs` take a single value.
#'      All other params are numeric vectors containing either one value or `niter` values.
#'      If one value is supplied, that value is used for all iterations of the model.
#'      If `niter`values are supplied, a different sequential value in the `niter` vector is used for each iteration.
#'      
#' @examples 
#' # Create a set of params
#' params <- roleParams(individuals_local = 100, individuals_meta = 1000,
#' species_meta = 10, speciation_local = 0.5, 
#' speciation_meta = 1, extinction_meta = 0.8, env_sigma = 0.5,
#' trait_sigma=1,comp_sigma = 0.1, dispersal_prob = 0.1, mutation_rate = 0.01,
#' equilib_escape = 1, num_basepairs = 250,
#' init_type = 'oceanic_island', niter = 2, niterTimestep = 2)
#' 
#' # Use it to create a roleModel
#' model <- roleModel(params)
#' @rdname roleParams
#' @export

setClass('roleParams',
         # slots that are 'functions' are allowed to time-vary 
         slots = c(
             # number of individuals and species in the local and meta
             individuals_local = 'function',
             individuals_meta = 'integer',
             species_meta = 'integer',
             
             # probs of speciation, extinction, dispersal
             speciation_local = 'function',
             speciation_meta = 'numeric',
             extinction_meta = 'numeric',
             dispersal_prob = 'function',
             
             env_optim = 'matrix',
             trait_sigma = 'numeric',
             env_sigma = 'numeric',
             comp_sigma = 'numeric',
             neut_delta = 'numeric',
             env_comp_delta = 'numeric',
             
             # gene evolution simulations
             mutation_rate = 'numeric',
             equilib_escape = 'numeric',
             alpha = 'function',
             num_basepairs = 'integer',
             
             init_type = 'character', 
             
             # iterations
             niter = 'integer', 
             niterTimestep = 'integer'
         )
)


#' @title Create a roleParams object 
#' @description Parameters include population sizes, rates of processes, the number of iterations
#' to run, and much more. 
#' 
#' @param individuals_local a single numeric or an "iter function"
#' @param individuals_meta a single integer
#' @param species_meta a single integer
#' 
#' @param speciation_local a single numeric or an "iter function"
#' @param speciation_meta a single numeric
#' @param extinction_meta a single numeric
#' @param dispersal_prob a single numeric or an "iter function"

#' @param trait_sigma a single numeric 
#' @param env_sigma a single numeric 
#' @param comp_sigma a single numeric 
#' @param neut_delta a single numeric 
#' @param env_comp_delta a single numeric 
#' 
#' @param mutation_rate a single numeric 
#' @param equilib_escape a single numeric 
#' @param alpha a single numeric or an "iter function"
#' @param num_basepairs a single integer
#' 
#' @param init_type a character vector either "oceanic_island" or "bridge_island"

#' @param niter an integer
#' @param niterTimestep an integer
#' 
#' @return a `roleParams` object 
#' 
#' @details `individuals_local`, `speciation_local` and `dispersal_prob`
#' are unique in that they are allowed to time-vary over the course of the model run. 
#' If you would like these to time-vary you can supply an "iter function" that takes a
#' model iteration number and returns the value that that param should take at that iteration in the model. 
#' 
#' @examples 
#' # create a default set of params but change the number of individuals in the local
#' p <- roleParams(individuals_local=500)
#' # create a new set of params but randomly deviating the speciation rate using an "iter function"
#' spfun <- function(i){rnorm(1,mean=0.1,sd=0.01)}
#' p <- roleParams(speciation_local=spfun)
#' @rdname roleParams
#' @export

roleParams <- function(individuals_local=100,
                       individuals_meta=1000,
                       species_meta=10,
                       
                       speciation_local=0,
                       speciation_meta=1,
                       extinction_meta=0.8,
                       dispersal_prob=0.01,
                       
                       env_optim = 0,
                       trait_sigma=1,
                       env_sigma=0,
                       comp_sigma=0,
                       neut_delta=1,
                       env_comp_delta=0.5,
                       
                       mutation_rate=0,
                       equilib_escape=0,
                       alpha=1,
                       num_basepairs=500,
                       
                       init_type='oceanic_island', 
                       niter=10, 
                       niterTimestep=NULL) {
    
    # if niterTimestep unspecified calculate one as rounded 1/10 of the iter 
    # plus 1
    if(is.null(niterTimestep)){
        niterTimestep <- as.integer(niter/10)
        if(niterTimestep <= 1){
            niterTimestep <- 2
        }
    }
    
    # check that iters and timesteps are correct
    if(niterTimestep > niter) {
        stop('`niter` cannot be less than `niterTimestep`')
    }
    
    # get a list of all the user supplied parameters 
    all_params <- as.list(environment())
    
    # get the types (i.e. 'function', 'numeric') of the slots of roleParam
    slot_types <- getSlots("roleParams")
    
    # get the names of each slot
    slot_names <- slotNames("roleParams")
    
    # deal with special case of `env_optim`
    all_params$env_optim <- matrix(all_params$env_optim, nrow = 1)
    
    # initialize new params object
    out_params <-  new("roleParams")
    
    # fill in the params object slots
    for(i in 1:length(all_params)) {
        # if slot type should be function, and function not provided, 
        # turn the user input into a function
        if(slot_types[i] == "function" & typeof(all_params[[i]]) != "closure") {
            if(length(all_params[[i]]) > 1) {
                stop("must provide a single value to `", names(all_params[i]), 
                     "`, or a function")
            } else {
                all_params[[i]] <- buildFun(all_params[[i]])
            }
        }
        
        val <- all_params[[i]] 
        
        # if slot needs an integer, coerce
        if(slot_types[i] == "integer") { 
            val <- as.integer(val)
        }
        
        # now value can be assigned to slot
        slot(out_params, names(all_params[i])) <- val 
    }
    
    return(out_params)
}


#' @title Wrapper around roleParams to create a "untb-flavored" (Unified Neutral Theory of Biodiversity) roleModel
#' @description Only arguments relevant to a UNTB neutral model is included
#' 
#' @param individuals_local  individuals_local
#'
#' @param individuals_meta individuals_meta
#' @param species_meta species_meta
#' @param speciation speciation
#' @param dispersal_prob dispersal_prob
#' @param init_type init_type
#' @param niter niter
#' @param niterTimestep niterTimestep
#' @return a `roleParams` object
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
        env_comp_delta = 1,
        dispersal_prob = dispersal_prob,
        mutation_rate = 1e-7,
        equilib_escape = 1,
        alpha = 10,
        num_basepairs = 250,
        init_type = init_type, 
        niter = niter,
        niterTimestep = niterTimestep))
}



# print method for `roleExperiment`
setMethod('show', signature = signature(object = 'roleParams'),
          definition = function(object) {
              n <- slotNames(object)
              tab <- data.frame(param = n, type = "function", 
                                starting_val = 1, ending_val = 1)
              
              for(i in 1:length(n)) {
                  s <- slot(object, n[i])
                  
                  if(is.function(s)) {
                      vals <- s(1:object@niter)
                      if(all(diff(vals) == 0)) {
                          tab[i, 2] <- "constant"
                          tab[i, 3:4] <- vals[1]
                      } else {
                          tab[i, 3:4] <- vals[c(1, object@niter)]
                      }
                  } else {
                      tab[i, 2] <- "constant"
                      tab[i, 3:4] <- s
                  }
              }
              
              print(tab, row.names = FALSE)
          }
)

# buildFun
# helper that, given a single value, builds a function that returns that value 
# stretched to niter in a vectorized fashion
# @param p something magical
#
# @return p stretched to niter

buildFun <- function(p) {
    p # what in the name of god why does this work
    f <- function(i) {
        out <- rep(p, length(i))
        
        return(out)
    }
    return(f)
}
