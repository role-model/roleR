#' @title One run of the RoLE model
#'
#' @description An S4 class to hold one model run of the RoLE model
#' 
#' @slot modelSteps a list of roleData objects, one for each recorded time step
#' @slot params a roleParams object with the RoLE model params
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
    
    
    locs <- localComm(indSpecies=rep(1,J),
                      indTrait = seq(1, 1.2, length.out = J),
                      indSeqs = rep('ATCG', J),
                      spGenDiv = c(1))
                      
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
    
    meta <- metaComm(spAbund =(Sm:1) / Sm, spTrait = 1:Sm)
    
    phylo <- ape::rphylo(Sm, params@speciation_meta, params@extinction_meta)
    phylo <- as(phylo, 'rolePhylo')
    
    dat <- roleData(localComm = locs, metaComm = meta, phylo = phylo)
    
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

