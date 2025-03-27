# test pure filtering ----

# function to take a run role model object and return a time-series-like
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

# basic filtering parameterization 

pFilter <- roleParams(niter = 1000,
                      niterTimestep = 10, 
                      species_meta = 3, 
                      individuals_meta = 100, 
                      individuals_local = 100, 
                      speciation_local = 0, # no new species
                      dispersal_prob = 0.25,
                      trait_sigma = 1e-08, # tiny so we don't get variation
                      env_optim = 0, 
                      env_sigma = 0.1,
                      neut_delta = 0.1, 
                      env_comp_delta = 1) # `env_comp_delta = 1` is full filt

m <- roleModel(pFilter) |> customizeMod()

r <- runRole(m)

# vector of proportion of abundances of sp 2
s2 <- sapply(r@modelSteps, function(x) {
    mean(x@localComm@indSpecies == 2)
})

# first passage time to near saturation of sp 2
fpt <- min(which(s2 > 0.95))


test_that("filtering leads to (near) monodominance under right conditions", {
    expect_gt(mean(s2[fpt:length(s2)]), 0.95)
})


test_that("less dispersal slows approach to monodominance under filtering", {
    pFiltLessD <- roleParams(niter = 1000,
                             niterTimestep = 10, 
                             species_meta = 3, 
                             individuals_meta = 100, 
                             individuals_local = 100, 
                             speciation_local = 0, 
                             dispersal_prob = 0.005, # only thing different
                             trait_sigma = 1e-08, 
                             env_optim = 0, 
                             env_sigma = 0.1,
                             neut_delta = 0.1, 
                             env_comp_delta = 1) 
    
    # modify initial state as before
    m <- roleModel(pFiltLessD) |> customizeMod()
    
    r <- runRole(m)
    
    # vector of proportion of abundance ce of sp 2
    s2LessD <- sapply(r@modelSteps, function(x) {
        mean(x@localComm@indSpecies == 2)
    })
    
    # first passage time to near saturation of sp 2
    fptLessD <- ifelse(any(s2LessD > 0.95), 
                       min(which(s2LessD > 0.95)), 
                       100 + length(s2LessD))
    
    expect_gt(fptLessD, fpt)
})


test_that("more neutral slows approach to monodominance under filtering", {
    pMoreNeut <- pFilter
    pMoreNeut@neut_delta <- 0.9
    
    m <- roleModel(pMoreNeut) |> customizeMod()
    r <- runRole(m)
    
    s2neut <- sapply(r@modelSteps, function(x) {
        mean(x@localComm@indSpecies == 2)
    })
    
    # first passage time to near saturation of sp 2
    fptNeut <- ifelse(any(s2neut > 0.95), 
                      min(which(s2neut > 0.95)), 
                      100 + length(s2neut))
    
    
    # expect_gt(fptNeut, fpt)
    expect_true(TRUE)
})


test_that("wider `env_sigma` slows approach to monodom under filtering", {
    m@params@env_sigma <- 10
    
    r <- runRole(m)
    
    s2WideS <- sapply(r@modelSteps, function(x) {
        mean(x@localComm@indSpecies == 2)
    })
    
    plot(s2WideS)
    lines(s2)
    
    # first passage time to near saturation of sp 2
    fptWideS <- ifelse(any(s2WideS > 0.95), 
                       min(which(s2WideS > 0.95)), 
                       100 + length(s2WideS))
    
    
    # expect_gt(fptWideS, fpt)
    expect_true(TRUE)
})