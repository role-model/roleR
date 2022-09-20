#' @title One run of the RoLE model
#'
#' @description An S4 class to hold one model run of the RoLE model
#' 
#' @details A model is first initialized using a set of params, then iterated using iterModel(modeL)
#' roleExperiments consist of many roleModels with different parameters
#' 
#' @slot modelSteps a list of roleData objects, one for each recorded time step
#' `niterTimestep` param defines recording interval
#' @slot params a `roleParams` object defining the model params
#' 
#' @examples 
#' Create a model using a default set of params, then run it
#' model <- roleModel(roleParams())
#' run <- iterModel(model)
#' 
#' @rdname roleModel
#' @export
#' 
setClass('roleModel',
         slots = c(params = 'roleParams', modelSteps = 'list'))

# constructor for roleModel
#' @rdname roleModel
#' @export

roleModel <- function(params) {
    J <- params@individuals_local[1]
    # if(J < 100) {
    #     stop('`individuals_local` (set in `roleParams`) cannot be less than 100')
    # }
    # 
    Sm <- params@species_meta
    # if(Sm < 200) {
    #     stop('`species_meta` (set in `roleParams`) cannot be less than 200')
    # }
    
    

                      
    #indSpTrt <- matrix(c(rep(1, J), seq(1, 1.2, length.out = J)), ncol = 2)
    
    #locs <- localComm(indSppTrt = indSpTrt, 
    #                  indSeqs = matrix(rep('ATCG', J), ncol = 1), 
    #                  sppGenDiv = matrix(1, ncol = 1))
    
    # for(s in 1:length(species)){
    #     sp <- species[1];
    #     abund <- 0
    #     traits <- c()
    #     for(r in 1:nrow(indSpTrt)){
    #         print(indSpTrt[1,1])
    #         print(sp)
    #         if(indSpTrt[r,1] == sp){
    #             abund <- abund + 1
    #             traits <- c(traits,indSpTrt(r,1))
    #         }
    #     }
    #     spAbundTrt(sp,0) <- abund;
    #     spAbundTrt(sp,1) <- mean(traits); 
    # }
    # locs@spAbundTrt <- spAbundTrt
    
    
    
    phylo <- ape::rphylo(Sm, params@speciation_meta, params@extinction_meta)
    
    meta <- metaComm(spAbund = .lseriesFromSN(params@species_meta, 
                                              params@individuals_meta), 
                     spTrait = ape::rTraitCont(phylo, 
                                               sigma = params@trait_sigma))
    
    # initialize indSpecies from random draw from meta (based on oceanic or
    # bridge island)
    if(params@init_type == 'oceanic_island') {
        initSpp <- rep(sample(params@species_meta, 1, 
                              prob = meta@spAbund), 
                       J)
    } else if(params@init_type == 'bridge_island') {
        initSpp <- sample(params@speciation_meta, J, 
                          replace = TRUE, prob = meta@spAbund)
    } else {
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
    if(niter > 100) {
        #stop('`niter` (set in `roleParams`) cannot be greater than 100')
    }
    
    niterTimestep <- params@niterTimestep
    
    # output data
    modelSteps <- vector('list', length = niter / niterTimestep + 1)
    modelSteps[[1]] <- dat
    #j <- 2 # counter to keep track of where to save data
    
    #for(i in 1:niter) {
        # update local comm
    #    locs@indSppTrt[i, ] <- c(i + 1, meta@sppAbundTrt[i + 1, 2])
    #    locs@sppGenDiv <- rbind(locs@sppGenDiv, 
    #                            matrix(1 / (i + 1), nrow = 1, ncol = 1))
        
        # write data every `niterTimestep`
    #    if(i %% niterTimestep == 0) {
            # over-write local comm in dat
    #        dat@localComm <- locs
            
            # save it
    #        modelSteps[[j]] <- dat
    #        j <- j + 1
    #    }
    #}
    
    # save last step if haven't already
    #if(i %% niterTimestep != 0) {
    #    modelSteps <- c(modelSteps, dat)
    #}
    
    return(new('roleModel', 
               params =  params, 
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
    thisSAD <- pika::sad(model = 'fish', par = beta)
    thisSAD <- pika::sad2Rank(thisSAD, S = S)
    
    # return relative abundances
    return(thisSAD / sum(thisSAD))
}


