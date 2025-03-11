# make sure traits are updated correctly (e.g. new trait is closer to parent 
# trait than other traits)

# test trait update with immigration ----

n <- 1
indsLocs <- 2
p <- roleParams(niter = n, 
                niterTimestep = 1,
                species_meta = 2, 
                individuals_meta = 100, 
                individuals_local = indsLocs, 
                dispersal_prob = 1,
                speciation_local = 0,
                comp_sigma = 0.5,
                env_sigma = 0.5,
                neut_delta = 1, # all neutral
                env_comp_delta = 1)

m <- roleModel(p)

# update the model to specific config for testing:
# in metacomm, only sp 1 is present
m@modelSteps[[1]]@metaComm@spAbund <- c(1, 0)

# in metacomm, spp have very different traits
metaTrt <- matrix(c(-100, 100), ncol = 1)
m@modelSteps[[1]]@metaComm@spTrait <- metaTrt

# init of local comm is all sp 2
m@modelSteps[[1]]@localComm@indSpecies <- c(2, 2)
m@modelSteps[[1]]@localComm@indTrait <- matrix(c(100, 100), ncol = 1)

# now we can test if trait of new individual (which *has* to be of sp 1) is 
# closer to trait of sp 1, as it should be
mr <- runRole(m)

IDs <- mr@modelSteps[[2]]@localComm@indSpecies
trts <- mr@modelSteps[[2]]@localComm@indTrait

# the test:
abs(trts[IDs == 1] - metaTrt[1, 1]) < abs(trts[IDs == 1] - metaTrt[2, 1])




# test trait update with birth ----

n <- 1
indsLocs <- 2
p <- roleParams(niter = n, 
                niterTimestep = 1,
                species_meta = 2, 
                individuals_meta = 100, 
                individuals_local = indsLocs, 
                dispersal_prob = 0,
                speciation_local = 0,
                comp_sigma = 0.5,
                env_sigma = 0.5,
                neut_delta = 1, # all neutral
                env_comp_delta = 1)

m <- roleModel(p)

# update the model to specific config for testing:

# spp have very different traits
metaTrt <- matrix(c(-100, 100), ncol = 1)

# init of local comm is one of each spp
m@modelSteps[[1]]@localComm@indSpecies <- c(1, 2)
m@modelSteps[[1]]@localComm@indTrait <- matrix(c(-100, 100), ncol = 1)

# now we can test if trait of new individual is closer to its parent 
mr <- runRole(m)

trtsInit <- mr@modelSteps[[1]]@localComm@indTrait
trtsFinal <- mr@modelSteps[[2]]@localComm@indTrait

newID <- which(trtsInit[, 1] != trtsFinal[, 1])
parentSp <- mr@modelSteps[[2]]@localComm@indSpecies[newID]
oldSp <- mr@modelSteps[[1]]@localComm@indSpecies[newID]


# the test:
abs(metaTrt[parentSp, 1] - trtsFinal[newID, 1]) < # diff of parent and offspring
    abs(diff(metaTrt[, 1])) # diff of two spp
