# test pure filtering ----

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

test_that("basic filtering model runs without error", {
    expect_no_error(runRole(roleModel(pFilter)))
})


# save this for further tests
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
    pMoreNeut@neut_delta <- 0.5
    
    m <- roleModel(pMoreNeut) |> customizeMod()
    r <- runRole(m)
    
    s2neut <- sapply(r@modelSteps, function(x) {
        mean(x@localComm@indSpecies == 2)
    })
    
    # first passage time to near saturation of sp 2
    fptNeut <- ifelse(any(s2neut > 0.95), 
                      min(which(s2neut > 0.95)), 
                      100 + length(s2neut))
    
    expect_gt(fptNeut, fpt)
})


test_that("when there is trait evolution, populations move toward envOpt", {
    pTrtEvo <- roleParams(niter = 100000,
                          niterTimestep = 1000, 
                          species_meta = 3, 
                          individuals_meta = 100, 
                          individuals_local = 100, 
                          speciation_local = 0, # no new species
                          dispersal_prob = 0.01, # immigration erases trt evo
                          trait_sigma = 50, # big
                          env_optim = 0, 
                          env_sigma = 0.1,
                          neut_delta = 0.1, 
                          env_comp_delta = 1)
    
    mTrtEvo <- roleModel(pTrtEvo) |> customizeMod()
    rTrtEvo <- runRole(mTrtEvo)
    
    trtEvo <- sapply(rTrtEvo@modelSteps, function(x) {
        spp <- factor(x@localComm@indSpecies, 1:x@phylo@n)
        tt <- tapply(x@localComm@indTrait, spp, mean)
        tt
    })

    # gather info so we can compare ancestral trait and new trait to optim
    ancTrt <- mTrtEvo@modelSteps[[1]]@metaComm@spTrait[1, 1]
    newTrt <- mean(trtEvo[1, ], na.rm = TRUE)
    e <- 0 # hard coding the env optim, might cause breaking
    
    expect_gt(abs(ancTrt - e), abs(newTrt - e))
})


test_that("no error occurs when one sp have trait exactly = optim", {
    # update basic filtering model so sp 2 is exactly at optim
    m@modelSteps[[1]]@metaComm@spTrait[2, 1] <- 0
    
    expect_no_error(runRole(m))
})

## behavior of wider sig is not what I thought, so leave this one alone for now
# test_that("wider `env_sigma` slows approach to monodom under filtering", {
#     pWider <- pFilter
#     pWider@env_sigma <- 10
# 
#     m <- roleModel(pWider) |> customizeMod()
#     r <- runRole(m)
#     
#     s2WideS <- sapply(r@modelSteps, function(x) {
#         mean(x@localComm@indSpecies == 2)
#     })
#     
#     # first passage time to near saturation of sp 2
#     fptWideS <- ifelse(any(s2WideS > 0.95), 
#                        min(which(s2WideS > 0.95)), 
#                        100 + length(s2WideS))
#     
#     
#     expect_gt(fptWideS, fpt)
# })

