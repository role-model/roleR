#' @title Parameters of one roleModel
#' @description An S4 class containing params for population sizes, rates of processes, the number of iterations
#' to run, and much more 
#' 
#' @slot individuals_local number of individuals in local community (J)
#' Determines the length of individual-level data vectors
#' When the model is initialized, all individuals are of a single species from the metacommunity
#' @slot individuals_meta number of individuals in meta community
#' Used in generating a log series of the initial abundances of species in the meta
#' @slot species_meta number of species in meta community
#' Determines the initial size of the phylogeny
#' @slot speciation_local probability of speciation occurring in the local comm each time step
#' The new local species can be either from a birth in the local comm or an immigration from the meta comm
#' @slot speciation_meta rate of speciation in meta community
#' This is used during the simulation of the start phylogeny before the model is run
#' The sum of speciation_meta and extinction_meta is the average lifetime of phylo branches, and the larger this value the less new individual traits will deviate
#' @slot extinction_meta rate of extinction in meta community
#' Like speciation_meta, used in the starting phylo simulation and in relation to traits
#' @slot trait_sigma rate of Brownian trait evolution in the meta community
#' Determines how much the trait of a new individual deviates from its parent; how fast traits change
#' @slot env_sigma selectivity of environmental filter; how strongly the environment selects which trait values are a good match for it
#' The larger the value, the less chance an individual will survive if it's far from the trait optimum (which is 0)
#' @slot comp_sigma selectivity of competition
#' @slot dispersal_prob probability of dispersal (immigration) occurring from the meta to the local
#' Every time step, either birth or immigration happens, so the probability of birth is 1 minus the dispersal_prob
#' 
#' @slot mutation_rate rate of sequence mutation to use in genetic simulations
#' @slot equilib_escape proportion of equilibrium required to halt the model as it is running and return it
#' @slot num_basepairs number of basepairs to use in genetic simulations
#' 
#' @slot init_type the biological model used to initialize; a single character string that can be either "oceanic_island", "bridge_island", or "bare_island"
#' The bridge island model has the initial individuals in the local comm arriving through a land bridge, while the oceanic has no bridge and is populated by a single dispersal
#' Thus in oceanic island all individuals are of a SINGLE species sampled proportional to meta comm species abunds, 
#' while in bridge island species individuals are sampled of MANY species proportional to their abundances
#' Bare island is related to oceanic, but instead of starting with "individuals_local" individuals of the sole sampled species, only 1 individual 
#' of that species appears and the rest of the space is filled with placeholder "rocks" representing unfilled space
#' @slot niter an integer specifying the number of time steps for the model to run
#' @slot niterTimestep an integer specifying the frequency (in numbers of 
#'     iterations) at which the model state is snapshotted and saved in a model's model steps object
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
                        # slots that are 'functions' are allowed to time-vary across the simulation - all others are not
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


#' @title Create a roleParams object with 
#' @description 
#' @param individuals_local a single numeric value or an iter function
#' @param individuals_meta number of individuals in meta community
#' Used in generating a log series of the initial abundances of species in the meta
#' @param species_meta number of species in meta community
#' Determines the initial size of the phylogeny
#' @param speciation_local probability of speciation occurring in the local comm each time step
#' The new local species can be either from a birth in the local comm or an immigration from the meta comm
#' @param speciation_meta rate of speciation in meta community
#' This is used during the simulation of the start phylogeny before the model is run
#' The sum of speciation_meta and extinction_meta is the average lifetime of phylo branches, and the larger this value the less new individual traits will deviate
#' @param extinction_meta rate of extinction in meta community
#' Like speciation_meta, used in the starting phylo simulation and in relation to traits
#' @param trait_sigma rate of Brownian trait evolution in the meta community
#' Determines how much the trait of a new individual deviates from its parent; how fast traits change
#' @param env_sigma selectivity of environmental filter; how strongly the environment selects which trait values are a good match for it
#' The larger the value, the less chance an individual will survive if it's far from the trait optimum (which is 0)
#' @param comp_sigma selectivity of competition
#' @param dispersal_prob probability of dispersal (immigration) occurring from the meta to the local
#' Every time step, either birth or immigration happens, so the probability of birth is 1 minus the dispersal_prob
#' 
#' @param mutation_rate rate of sequence mutation to use in genetic simulations
#' @param equilib_escape proportion of equilibrium required to halt the model as it is running and return it
#' @param num_basepairs number of basepairs to use in genetic simulations
#' 
#' @param init_type the biological model used to initialize; a single character string that can be either "oceanic_island" or "bridge_island"
#' The bridge island model has the initial individuals in the local comm arriving through a land bridge, while the oceanic has no bridge and is populated by a single dispersal
#' Thus in oceanic island all individuals are of a SINGLE species sampled proportional to meta comm species abunds, 
#' while in bridge island species individuals are sampled of MANY species proportional to their abundances
#' @param niter an integer specifying the number of time steps for the model to run
#' @param niterTimestep an integer specifying the frequency (in numbers of 
#'     iterations) at which the model state is snapshotted and saved in a model's model steps object
#' @return a `roleParams` object 
#' 
#' @rdname roleParams
#' @export

roleParams <- function(individuals_local=100,
                       individuals_meta=1000,
                       species_meta=10,
                       
                       speciation_local=0,
                       speciation_meta=1,
                       extinction_meta=0.8,
                       dispersal_prob=0.1,
                       
                       trait_sigma=1,
                       env_sigma=0,
                       comp_sigma=0,
                       neut_delta=1,
                       env_comp_delta=0.5,
                       
                       mutation_rate=0,
                       equilib_escape=0,
                       alpha=0,
                       num_basepairs=0,
                       
                       init_type='oceanic_island', 
                       niter=10, 
                       niterTimestep=NULL) {
    
    # if niterTimestep unspecified calculate one as rounded 1/10 of the iter plus 1
    if(is.null(niterTimestep)){
        niterTimestep <- as.integer(niter/10)
        if(niterTimestep <= 1){
            niterTimestep <- 2
        }
    }
    
    # check that iters and timesteps are correct
    if(!niter%%1==0 | !niterTimestep%%1==0){ # check integer
        stop('niter and niterTimestep must be numeric integers (cannot be decimal)')
    }
    if(length(niter) > 1 | length(niterTimestep) > 1 | niterTimestep > niter) {
        stop('must supply a single value for `niter`and niterTimestep, and niter cannot be less than niterTimestep')
    }
    
    # get a list of all the user supplied parameters 
    # old way of doing this is all_params <- list(individuals_local,...
    # all_params <- as.list(environment())
    all_params <- list(individuals_local,
                       individuals_meta,
                       species_meta,
                       
                       speciation_local,
                       speciation_meta,
                       extinction_meta,
                       dispersal_prob,
                       
                       trait_sigma,
                       env_sigma,
                       comp_sigma,
                       neut_delta,
                       env_comp_delta,
                       
                       mutation_rate,
                       equilib_escape,
                       alpha,
                       num_basepairs,
                       
                       init_type, 
                       niter, 
                       niterTimestep)
    
    # get the types (i,e. 'function','numeric') of the slots of the roleParams class
    slot_types <- getSlots("roleParams")
    # get the names of each slot
    slot_names <- slotNames("roleParams")
    
    
    # for every slot in the list of slots
    for(i in 1:length(all_params)){

        # if the slot type is a function, and the user input is NOT a function...
        if(slot_types[i] == "function" & (typeof(all_params[[i]]) != "closure")){
            # replace the single user-supplied value with the function
            all_params[[i]] <- buildFun(all_params[[i]])
        }
    }
    
    # singleValParams <- c('individuals_meta', 'species_meta',
    #                      'speciation_meta', 'extinction_meta', 'trait_sigma', 'env_sigma', 'comp_sigma',
    #                      'equilib_escape', 'num_basepairs', 'init_type', 
    #                      'niter', 'niterTimestep', 'neut_delta', 'env_comp_delta')
    
    # create params to return and populate with updated values (values that are replaced with functions)
    out_params <-  new('roleParams')
    for(i in 1:length(getSlots('roleParams'))){ # for each slot
        val <- all_params[[i]] # get value to assign
        if(slot_types[i] == "integer"){ # if slot needs an integer, coerce
            val <- as.integer(val)
        }
        slot(out_params,slot_names[i]) <- val # add the value to the corresponding slot name in out
    }
    
    return(out_params)
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
        env_comp_delta = 1,
        dispersal_prob = dispersal_prob,
        mutation_rate = 0.01,
        equilib_escape = 1,
        alpha = 10,
        num_basepairs = 250,
        init_type = init_type, 
        niter = niter,
        niterTimestep = niterTimestep))
}

# constructor
#' @rdname untbParams
#' @export

sp1_gr = 0.1
sp2_gr = 0.6
sp1_k = 50
sp2_k = 70
alpha12 = 0.1
alpha21 = 0.2
niter = 1000
lvParams <- function(sp1_n,sp2_n,sp1_gr,sp2_gr,sp1_k,sp2_k,alpha12,alpha21,niter) {
    max_gr_sp <- which.max(c(sp1_gr,sp2_gr))
    max_gr <- max(sp1_gr,sp2_gr)
    
    mu10 <- ifelse(max_gr_sp == 1, 0, max_gr - sp1_gr)
    mu20 <- ifelse(max_gr_sp == 2, 0, max_gr - sp2_gr)
    
    la10 <- sp1_gr + mu10
    la20 <- sp2_gr + mu20
    
    mu11 <- sp1_gr / sp1_k
    mu22 <- sp2_gr / sp2_k
    
    mu12 <- sp1_gr * alpha12 / sp1_k
    mu21 <- sp2_gr * alpha21 / sp2_k
    
    return(roleParams(
        individuals_local = sp1_k + sp2_k, # J (number of indv + rocks) = the total carrying capacity
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
        mutation_rate = 0.01,
        equilib_escape = 1,
        alpha = 10,
        num_basepairs = 250,
        init_type = init_type, 
        niter = niter,
        niterTimestep = niterTimestep))
}

# helper that, given a single value, builds a function
#   that returns that value stretched to niter in a vectorized fashion
buildFun <- function(p) {
    p # what in the name of god why does this work
    f <- function(i) {
        out <- rep(p, length(i))
        
        return(out)
    }
    return(f)
}