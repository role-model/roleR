#' @title A single RoLE model.
#'
#' @description An S4 class that holds one RoLE model realization. A model is 
#'     first initialized using a set of parameters, then run using those 
#'     parameters.
#' 
#' @slot modelSteps A list of `roleData` objects, one for each saved snapshot.
#'     For example, the 3rd saved snapshot is accessed at `modelSteps[[3]]`.
#'     Models that are not yet run only have one timestep, i.e. `modelSteps[[1]]`
#' @slot params A `roleParams` object containing the params of the model
#' @slot info a `data.frame` with one row for each saved snapshot (unrun models
#'     will only have 1 row); columns are parameters, cells are parameter values 
#'     at each snapshot
#' 
#' @details See the `roleR_intro` vignette for an example modeling workflow.
#' @examples 
#' # Create a model using a default set of params, then run it.
#' m <- roleModel(roleParams())
#' m <- runRole(m)
#' 
#' @rdname roleModel-class
#' @export

setClass('roleModel',
         slots = c(params = 'roleParams', 
                   info = 'data.frame',
                   modelSteps = 'list'))

#' @title Create a roleModel.
#' @param params The params, an object of class `roleParams`, to use when the 
#'     model is run.
#' @return A ready-to-run `roleModel`.
#' 
#' @rdname roleModel
#' @export

roleModel <- function(params) {
    J <- params@individuals_local(1)
    Sm <- params@species_meta

    phylo <- ape::rphylo(Sm, params@speciation_meta, params@extinction_meta)

    meta <- metaComm(spAbund = .lseriesFromSN(params@species_meta, 
                                              params@individuals_meta), 
                     spTrait = ape::rTraitCont(phylo, 
                                               sigma = params@trait_sigma))
    
    # initialize indSpecies from random draw from meta (based on oceanic or
    # bridge island)
    if(params@init_type == 'oceanic_island'){
        initSpp <- rep(sample(params@species_meta, 1, 
                              prob = meta@spAbund), 
                       J)
    } 
    else if(params@init_type == 'bridge_island'){
        initSpp <- sample(params@species_meta, J, 
                          replace = TRUE, prob = meta@spAbund)
    } 
    else if(params@init_type == 'bare_island'){
        initSpp <- rep(sample(params@species_meta, 1, 
                              prob = meta@spAbund), 
                       J)
        initSpp[2:length(initSpp)] <- 0
    }
    else{
        stop('`init_type` must be one of "oceanic_island" or "bridge_island"')
    }

    # initialize traits based on spp ID
    initTrait <- meta@spTrait[initSpp]
    
    locs <- localComm(indSpecies = initSpp,
                      indTrait = initTrait,
                      indSeqs = character(J),
                      spGenDiv = numeric(0))
    
    dat <- roleData(localComm = locs, 
                    metaComm = meta, 
                    phylo = as(phylo, 'rolePhylo'))
    
    niter <- params@niter
    niterTimestep <- params@niterTimestep
    
    # calculate generations
    allTStep <- 1:(niter + 1)
    allJ <- params@individuals_local(allTStep)
    gens <- (allTStep - 1) * 2 / allJ # scale so generations start at 0
    
    # output data
    modelSteps <- vector('list', length = niter / niterTimestep + 1)
    modelSteps[[1]] <- dat # initial state 
    
    
    
    
    # create info data.frame
    tstep <- seq(1, params@niter + 1, params@niterTimestep)
    info <- as.data.frame(getValuesFromParams(params, tstep))
    info <- info[, !grepl('niter', names(info))]
    info$timestep <- tstep - 1 # scale so initial state is tstep 0
    info$generations <- gens[tstep]
    
    
    return(new('roleModel', 
               params =  params, 
               info = info,
               modelSteps = modelSteps))
}


# ----
#' @description function to solve for parameter of logseries
#' @param S number of species
#' @param N number of individuals

.lseriesFromSN <- function(S, N) {
    # solve for alpha paramter
    # browser()
    asol <- uniroot(interval = c(.Machine$double.eps^0.25,
                                 .Machine$integer.max),
                    f = function(a) {
                        a * log(1 + N / a) - S
                    })
    
    # calculate p parameter and beta (as used by pika)
    p <- 1 - exp(-S / asol$root)
    beta <- -log(p)
    
    # calculate idealized SAD from parameter
    thisSAD <- pika::sad(model = 'lseries', par = beta)
    thisSAD <- pika::sad2Rank(thisSAD, S = S)
    
    # return relative abundances
    return(thisSAD / sum(thisSAD))
}

# setValidity()
