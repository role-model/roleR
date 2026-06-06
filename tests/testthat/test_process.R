# test competition modulated by traits ----

test_that("intraspecific trait variation increases abundance of that species", {
    # sp1: all indis at the same trait = heavy within-species competition
    # sp2: inds spread across trait space = low within-species competition

    n <- 5000
    J <- 100
    p <- roleParams(
        niter = n,
        niterTimestep = n,
        species_meta = 2,
        individuals_meta = 100,
        individuals_local = J,
        dispersal_prob = 0, # closed community: no immigration 
        speciation_local = 0,
        trait_sigma = 1e-8, # traits barely evolve, initial variation preserved
        comp_sigma = 0.5,
        env_sigma = 1,
        neut_delta = 0,
        env_comp_delta = 0 # pure competition
    )

    m <- roleModel(p)

    # equal starting abundances: 50 each
    m@modelSteps[[1]]@localComm@indSpecies <- c(rep(1, J / 2), rep(2, J / 2))
    m@modelSteps[[1]]@localComm@spAbund <- rep(J / 2, 2)

    m@modelSteps[[1]]@localComm@spAbundHarmMean <- rep(J / 2, 2)
    m@modelSteps[[1]]@localComm@spPropHarmMean <- rep(0.5, 2)
    m@modelSteps[[1]]@localComm@spLastOriginStep <- rep(0, 2)
    m@modelSteps[[1]]@localComm@spLastOriginGen <- rep(0, 2)
    
    # sp1: zero variation (all at trait = 0), maximum intra-sp comp
    # sp2: high variation (spread from -5 to +5), less intra-sp comp
    trt1 <- rep(0, J / 2)
    trt2 <- seq(-5, 5, length.out = J / 2)
    m@modelSteps[[1]]@localComm@indTrait <- matrix(c(trt1, trt2), 
                                                   ncol = 1)

    mr <- runRole(m)
    
    # sp2 should end up more abundant than sp1
    expect_gt(diff(mr@modelSteps[[2]]@localComm@spAbund), 0)
})


# test filtering (closest to optimum wins) ----

test_that("filtering: species closest to env optim takes over", {
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

    # spp 1 and 2 have equal abundance, sp 3 has 0 abund
    m@modelSteps[[1]]@metaComm@spAbund <- c(0.5, 0.5, 0)

    # sp 1 matches optim perfectly, spp 2 and 3 are far from optim
    m@modelSteps[[1]]@metaComm@spTrait <- matrix(c(0, 10, 10), ncol = 1)

    # init: all sp 3
    m@modelSteps[[1]]@localComm@indSpecies <- rep(3, indsLocs)
    m@modelSteps[[1]]@localComm@indTrait <- matrix(10, nrow = indsLocs, ncol = 1)

    mr <- runRole(m)

    # sp 1 should almost completely take over; allow 1 straggler (last event
    # might be a dispersal of sp 2 which won't have time to die)
    expect_gte(sum(mr@modelSteps[[3]]@localComm@indSpecies == 1), indsLocs - 1)
})


test_that("filtering: two spp with same trait have ~equal abundance", {
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
    m@modelSteps[[1]]@metaComm@spAbund <- c(0.5, 0.5, 0)

    # spp 1 and 2 both match optim; should have similar abundance
    m@modelSteps[[1]]@metaComm@spTrait <- matrix(c(0, 0, 10), ncol = 1)
    m@modelSteps[[1]]@localComm@indSpecies <- rep(3, indsLocs)
    m@modelSteps[[1]]@localComm@indTrait <- matrix(10, nrow = indsLocs, ncol = 1)

    mr <- runRole(m)
    nsp1 <- sum(mr@modelSteps[[3]]@localComm@indSpecies == 1)

    # test uses binomial assumption for distribution of sp1 count
    expect_lte(nsp1, qbinom(1 - .Machine$double.eps^0.5, 100, 0.5))
    expect_gte(nsp1, qbinom(1 - .Machine$double.eps^0.5, 100, 0.5, lower.tail = FALSE))
})



# test speciation happens only when it should ----

test_that("speciation increases species richness", {
    n <- 10000
    nspp0 <- 2
    p <- roleParams(niter = n,
                    niterTimestep = n / 2,
                    species_meta = nspp0,
                    individuals_meta = 100,
                    individuals_local = 100,
                    dispersal_prob = 0.25,
                    speciation_local = 0.001, # high speciation prob
                    comp_sigma = 0.5,
                    env_sigma = 0.5,
                    neut_delta = 0.5,
                    env_comp_delta = 1)

    m <- roleModel(p)
    
    mr <- runRole(m)

    # with high speciation rate, final nspp should exceed initial nspp
    expect_gt(mr@modelSteps[[length(mr@modelSteps)]]@phylo@n, 
              nspp0)
})




# test trait update with immigration ----

test_that("immigrating individual inherits parent species trait", {
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

    # only sp 1 in metacomm
    m@modelSteps[[1]]@metaComm@spAbund <- c(1, 0)

    # spp have very different traits
    metaTrt <- matrix(c(-100, 100), ncol = 1)
    m@modelSteps[[1]]@metaComm@spTrait <- metaTrt

    # init: all sp 2
    m@modelSteps[[1]]@localComm@indSpecies <- c(2, 2)
    m@modelSteps[[1]]@localComm@indTrait <- matrix(c(100, 100), ncol = 1)

    mr <- runRole(m)

    IDs <- mr@modelSteps[[2]]@localComm@indSpecies
    trts <- mr@modelSteps[[2]]@localComm@indTrait

    # the new sp 1 individual's trait should be closer to sp 1 trait than sp 2 trait
    expect_lt(abs(trts[IDs == 1] - metaTrt[1, 1]),
              abs(trts[IDs == 1] - metaTrt[2, 1])
    )
})


# test trait update with birth ----

test_that("offspring trait is closer to parent species trait than to other species", {
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

    # spp have very different traits
    metaTrt <- matrix(c(-100, 100), ncol = 1)

    # init: one of each spp
    m@modelSteps[[1]]@localComm@indSpecies <- c(1, 2)
    m@modelSteps[[1]]@localComm@indTrait <- matrix(c(-100, 100), ncol = 1)

    mr <- runRole(m)

    trtsInit <- mr@modelSteps[[1]]@localComm@indTrait
    trtsFinal <- mr@modelSteps[[2]]@localComm@indTrait

    newID <- which(trtsInit[, 1] != trtsFinal[, 1])
    parentSp <- mr@modelSteps[[2]]@localComm@indSpecies[newID]

    # offspring trait should be closer to its parent species trait than to the other species
    expect_lt(
        abs(metaTrt[parentSp, 1] - trtsFinal[newID, 1]),
        abs(diff(metaTrt[, 1]))
    )
})
