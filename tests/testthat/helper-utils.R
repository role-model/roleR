# a functio to take a run role model object and return a time-series-like
# matrix where rows are time steps and there is a column for each spp
# cells are proportional abundances

ats <- function(r) {
    nspp <- r@modelSteps[[length(r@modelSteps)]]@phylo@n
    
    res <- lapply(r@modelSteps, function(x) {
        sapply(1:nspp, function(i) {
            mean(x@localComm@indSpecies == i)
        })
    })
    
    return(do.call(rbind, res))
}

# function to customize model before running
# modify initial state so we have:
#    - one sp close to the optim
#    - equal abundances in metacomm
#    - initial sp is one not close to optim

customizeMod <- function(mod) {
    J <- mod@params@individuals_local(1)
    trts <- c(-3, 1e-06, 10)
    metaA <- rep(1/3, 3)
    initSp <- 1
    
    mod@modelSteps[[1]]@metaComm@spAbund <- metaA
    mod@modelSteps[[1]]@metaComm@spTrait <- matrix(trts)
    mod@modelSteps[[1]]@localComm@indSpecies <- rep(initSp, J) 
    mod@modelSteps[[1]]@localComm@indTrait <- matrix(trts[initSp], 
                                                     nrow = J, 
                                                     ncol = 1)
    
    return(mod)
}

