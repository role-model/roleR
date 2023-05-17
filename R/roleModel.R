
#' @title A single RoLE model.
#'
#' @description An S4 class that holds a RoLE eco-evolutionary process model.
#' A model is first initialized using a set of parameters, then run using those parameters.
#' 
#' @slot modelSteps A list of `roleData` objects, one for each snapshot of the model that were recorded as the model ran.
#' For example, the 3rd saved snapshot is accessed at modelSteps[[3]].
#' Models that are not yet run only have one timestep in modelSteps[[1]].
#' @slot params A `roleParams` object containing the params to use when the model is run.
#' 
#' @details See the `roleR_intro` vignette for an example modeling workflow.
#' @examples 
#' # Create a model using a default set of params, then run it.
#' m <- roleModel(roleParams())
#' m <- runRole(m)
#' @include roleModel.R roleParams.R
#' @rdname roleModel
#' @export
#' 
setClass('roleModel',
         slots = c(params = 'roleParams', modelSteps = 'list'))

#' @title Create a roleModel.
#' @param params The params to use when the model is run.
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
        stop('`init_type` must be one of `"oceanic_island"` or `"bridge_island"`')
    }
    
    # initialize traits based on spp ID
    initTrait <- meta@spTrait[initSpp]
    
    locs <- localComm(indSpecies = initSpp,
                      indTrait = initTrait,
                      indSeqs = rep('ATCG', J), # leave genetic stuff alone
                      spGenDiv = c(1))
    
    dat <- roleData(localComm = locs, 
                    metaComm = meta, 
                    phylo = as(phylo, 'rolePhylo'))
    
    niter <- params@niter
    niterTimestep <- params@niterTimestep
    
    # output data
    modelSteps <- vector('list', length = niter / niterTimestep + 1)
    modelSteps[[1]] <- dat
    
    return(new('roleModel', 
               params =  params, 
               modelSteps = modelSteps))
}



# non-exported helper function to solve for parameter of logseries
#' @param S number of species
#' @param N number of individuals

.lseriesFromSN <- function(S, N) {
    qfun <- function(p, beta) {
        extraDistr::qlgser(p, theta = exp(-beta))
    }
    
    # solve for alpha paramter
    asol <- uniroot(interval = c(.Machine$double.eps^0.25,
                                 .Machine$integer.max),
                    f = function(a) {
                        a * log(1 + N / a) - S
                    })
    
    # calculate p parameter and beta (as used by pika)
    p <- 1 - exp(-S / asol$root)
    beta <- -log(p)
    
    
    rank <- qfun(seq(1, 1/S, length = S) - 1/(2 * S), beta = beta)
    
    return(rank / sum(rank))
}


