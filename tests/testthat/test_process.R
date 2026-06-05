# test competition modulated by traits ----

test_that("intraspecific trait variation increases abundance of that species", {
    
})



# test filtering (closest to optimum wins) ----

n <- 10000
indsLocs <- 100
p <- roleParams(niter = n, 
                niterTimestep = n / 2,
                species_meta = 3, 
                individuals_meta = 100, 
                individuals_local = indsLocs, 
                dispersal_prob = 1,
                speciation_local = 0,
                comp_sigma = 0.5,
                env_sigma = 0.5,
                neut_delta = 0, 
                env_comp_delta = 1)

m <- roleModel(p)

# update the model to specific config for testing:
# in metacomm, spp 1 and 2 have equal abundance, sp 3 has 0 abund
m@modelSteps[[1]]@metaComm@spAbund <- c(0.5, 0.5, 0)

# in metacomm, spp 1 matches optim perfectly, spp 2 and 3 are far from optim
m@modelSteps[[1]]@metaComm@spTrait <- matrix(c(0, 10, 10), ncol = 1)

# init of local comm is all sp 3
m@modelSteps[[1]]@localComm@indSpecies <- rep(3, indsLocs)
m@modelSteps[[1]]@localComm@indTrait <- matrix(10, nrow = indsLocs, ncol = 1)

# now we can test if all of the sp 3 individuals are replaced by sp 1
# which is the sp that matches the optim
mr <- runRole(m)

# the test. total number of 1's can be `indsLocs - 1` because the last event
# might be a dispersal of sp 2 from the metacomm and as the last event, it 
# won't have time to die
sum(mr@modelSteps[[3]]@localComm@indSpecies == 1) >= indsLocs - 1

# now we can re-run the model but make sp 2 have same trait as sp 1 and we 
# should see them both have about equal abundance

# make spp 1 and 2 match
m@modelSteps[[1]]@metaComm@spTrait <- matrix(c(0, 0, 10), ncol = 1)

# re-run
mr <- runRole(m)

# test uses binomial assumption for dist of 1's and 2's
sum(mr@modelSteps[[3]]@localComm@indSpecies == 1) <= 
    qbinom(1 - .Machine$double.eps^0.5, 100, 0.5) &
    sum(mr@modelSteps[[3]]@localComm@indSpecies == 1) >= 
    qbinom(1 - .Machine$double.eps^0.5, 100, 0.5, lower.tail = FALSE)


# test that traits don't matter in neutrality ----
n <- 10000
indsLocs <- 100
p <- roleParams(niter = n, 
                niterTimestep = n / 2,
                species_meta = 2, 
                individuals_meta = 100, 
                individuals_local = indsLocs, 
                dispersal_prob = 1,
                speciation_local = 0,
                comp_sigma = 0.5,
                env_sigma = 0.5,
                neut_delta = 1, # ALL NEUTRAL DESPITE OTHER PARAMS 
                env_comp_delta = 1)

m <- roleModel(p)

# update the model to specific config for testing:
# in metacomm, spp 1 and 2 have equal abundance, sp 3 is nill
m@modelSteps[[1]]@metaComm@spAbund <- c(0.5, 0.5, 0)

# in metacomm, spp 1 matches optim perfectly, spp 2 and 3 are far from optim
m@modelSteps[[1]]@metaComm@spTrait <- matrix(c(0, 10, 10), ncol = 1)

# init of local comm is all sp 3
m@modelSteps[[1]]@localComm@indSpecies <- rep(3, indsLocs)
m@modelSteps[[1]]@localComm@indTrait <- matrix(10, nrow = indsLocs, ncol = 1)

# now we can test if sp 1 and sp 2 have approx same abund despite diff traits
mr <- runRole(m)

# the test. total number of 1's can be `indsLocs - 1` because the last event
# might be a dispersal of sp 2 from the metacomm and as the last event, it 
# won't have time to die
nsp1 <- sum(mr@modelSteps[[3]]@localComm@indSpecies == 1)
nsp2 <- sum(mr@modelSteps[[3]]@localComm@indSpecies == 2)

# test again uses binomial assumption and the fact that if:
# X ~ binom(n, 0.5) and Y ~ binom(n, 0.5), 
# then n - Y ~ binom(n, 0.5) and X - Y + n ~ binom(2 * n, 0.5)
nsp1 - nsp2 + indsLocs <= 
    qbinom(1 - .Machine$double.eps^0.5, indsLocs * 2, 0.5) &
    nsp1 - nsp2 + indsLocs >= 
    qbinom(1 - .Machine$double.eps^0.5, indsLocs * 2, 0.5, lower.tail = FALSE)



# test speciation happens only when it should ----

n <- 10000
indsLocs <- 100
p <- roleParams(niter = n, 
                niterTimestep = n / 2,
                species_meta = 2, 
                individuals_meta = 100, 
                individuals_local = indsLocs, 
                dispersal_prob = 0.25,
                speciation_local = 0.001, # high ass speciation prob for testing
                comp_sigma = 0.5,
                env_sigma = 0.5,
                neut_delta = 0.5, 
                env_comp_delta = 1)

m <- roleModel(p)
mr <- runRole(m)



# test speciation updates spp IDs correctly ----


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

