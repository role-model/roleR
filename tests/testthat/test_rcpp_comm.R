test_that("roleData and params extraction into rcpp works", {
    p <- roleParams(niter = 10, 
                    niterTimestep = 5,
                    species_meta = 5, 
                    individuals_meta = 100, 
                    individuals_local = 10, 
                    comp_sigma = 0.5,
                    neut_delta = 0.1, 
                    env_comp_delta = 0.5)
    
    m <- roleModel(p)
    
    expect_no_error(roleR:::roleCommTester(m@modelSteps[[1]], p))
})





# test pure competition ----

pComp <- roleParams(niter = 10000,
                    niterTimestep = 100, 
                    species_meta = 3, 
                    individuals_meta = 100, 
                    individuals_local = 100, 
                    speciation_local = 0, # no new species
                    dispersal_prob = 0.25,
                    trait_sigma = 1e-08, # tiny so we don't get variation
                    env_optim = 0, 
                    env_sigma = 0.1,
                    comp_sigma = 5,
                    neut_delta = 0.1, 
                    env_comp_delta = 0) # `env_comp_delta = 0` is full comp



mComp <- roleModel(pComp)

test_that("env optim doesn't matter under pure competition", {
    # modify initial state so we have:
    #    - one sp close to the optim (to check that optim no matter)
    #    - equal abundances in metacomm
    
    initSp <- 1
    trts <- c(0, -10, 10)
    
    mComp@modelSteps[[1]]@metaComm@spAbund <- rep(1/3, 3)
    mComp@modelSteps[[1]]@metaComm@spTrait <- matrix(trts)
    mComp@modelSteps[[1]]@localComm@indSpecies <- 
        rep(initSp, pComp@individuals_local(1)) 
    mComp@modelSteps[[1]]@localComm@indTrait <- 
        matrix(trts[initSp], nrow = pComp@individuals_local(1), ncol = 1)
    
    rComp <- runRole(mComp)
    
    # time steps 1--3 are burn in (based on ocular analysis)
    s1 <- sapply(rComp@modelSteps[-(1:3)], function(x) {
        mean(x@localComm@indSpecies == 1)
    })
    
    nrep <- length(s1)
    n <- rComp@params@individuals_local(1)
    p <- 1 / 3
    crit <- qnorm(0.9999, mean = p, sd = sqrt(p * (1 - p) / n / nrep)) - 1/3
    
    # we are testing if the mean proportion of sp 1 (which is the one closest
    # to the env optim) is substantially different from 1/3; 1/3 is the prop
    # it should be if env doesn't matter but competition does
    # we approximate the expected difference between the obs prop and 1/3
    # with a normal distribution
    expect_lt(abs(mean(s1) - 1/3), crit)
})

test_that("species with distinct trait is more abundant than others", {
    # modify initial state so we have:
    #    - one sp has distinct trait, others have same trait
    #    - equal abundances in metacomm
    
    initSp <- 1 # sp1 and sp2 have same trait, initialize with sp1
    trts <- c(-10, -10, 10)
    
    mComp@modelSteps[[1]]@metaComm@spAbund <- rep(1/3, 3)
    mComp@modelSteps[[1]]@metaComm@spTrait <- matrix(trts)
    mComp@modelSteps[[1]]@localComm@indSpecies <- 
        rep(initSp, pComp@individuals_local(1)) 
    mComp@modelSteps[[1]]@localComm@indTrait <- 
        matrix(trts[initSp], nrow = pComp@individuals_local(1), ncol = 1)
    
    rComp <- runRole(mComp)
    
    # time steps 1--3 are burn in (based on ocular analysis)
    a <- ats(rComp)[-(1:3), ]
    
    a <- apply(a, 2, function(x) {
        c(mean(x), mean(x) + sd(x) * c(-2, 2))
    })
    
    # test that abund of sp 1 and sp 2 are similar (w/n 2 sd)
    # and that abund of sp 3 is higher (above 2 sd)
    test1 <- a[1, 1] <= a[3, 2] & a[1, 1] >= a[2, 2]
    test2 <- a[1, 2] <= a[3, 1] & a[1, 2] >= a[2, 1]
    test3 <- mean(a[1, 1:2]) <= a[2, 3]
    
    expect_true(test1 & test2 & test3)
})




test_that("decreased width of comp kernal lessens impact of trait diffs", {
    # modify initial state so we have:
    #    - one sp has distinct trait, others have same trait
    #    - equal abundances in metacomm
    
    initSp <- 1 # sp1 and sp2 have similar trait, initialize with sp1
    trts <- c(-8, -10, 10)
    
    mComp@modelSteps[[1]]@metaComm@spAbund <- rep(1/3, 3)
    mComp@modelSteps[[1]]@metaComm@spTrait <- matrix(trts)
    mComp@modelSteps[[1]]@localComm@indSpecies <- 
        rep(initSp, pComp@individuals_local(1)) 
    mComp@modelSteps[[1]]@localComm@indTrait <- 
        matrix(trts[initSp], nrow = pComp@individuals_local(1), ncol = 1)
    
    # make another model object with narrow comp kernal
    mCompNrw <- mComp
    mCompNrw@params@comp_sigma <- 0.5
    
    rComp <- runRole(mComp)
    rCompNrw <- runRole(mCompNrw)
    
    # time steps 1--3 are burn in (based on ocular analysis)
    a <- ats(rComp)[-(1:3), ]
    
    a <- apply(a, 2, function(x) {
        c(mean(x), sd(x))
    })
    
    aNrw <- ats(rCompNrw)
    
    aNrw <- apply(aNrw, 2, function(x) {
        c(mean(x), sd(x))
    })
    
    # we expect the differences in abundances when comp kernal is narrow
    # to be less than difference sin abundances with wider
    aDiff <- (dist(a[1, ]) / mean(a[2, ])) |> mean()
    aDiffNrw <- (dist(aNrw[1, ]) / mean(aNrw[2, ])) |> mean()
    expect_lt(aDiffNrw, aDiff)
})



# tests with just one time step to see exactly who dies ----



# test pure neutrality ----

pNeut <- roleParams(niter = 10000,
                    niterTimestep = 10, 
                    species_meta = 3, 
                    individuals_meta = 100, 
                    individuals_local = 100, 
                    speciation_local = 0, # no new species
                    dispersal_prob = 0.025,
                    trait_sigma = 1e-08, # tiny so we don't get variation
                    env_optim = 0, 
                    env_sigma = 0.1,
                    neut_delta = 1, 
                    env_comp_delta = 1)

m <- roleModel(pNeut)

# modify initial state so we have:
#    - one sp close to the optim
#    - equal abundances in metacomm 
#    - initial sp is one not close to optim

m@modelSteps[[1]]@metaComm@spAbund <- rep(1/3, 3) 
m@modelSteps[[1]]@metaComm@spTrait <- matrix(c(-10, 0, -10))
m@modelSteps[[1]]@localComm@indSpecies <- rep(1, 100) # make sure this is sp 1
m@modelSteps[[1]]@localComm@indTrait <- matrix(-3, nrow = 100, ncol = 1)

rNeut <- runRole(m)


test_that("variance of abundances under neutral is greater than under comp", {
    aNeut <- ats(rNeut)
    aComp <- mComp |> runRole() |> ats()
    
    vNeut <- apply(aNeut, 2, var) |> mean()
    vComp <- apply(aComp, 2, var) |> mean()
    
    expect_gt(vNeut, vComp)
})


test_that("species closest to env optim is not dominant", {
    pNeut@niter <- 100000L
    pNeut@niterTimestep <- 500L
    pNeut@env_comp_delta <- 1 # as if this was filtering
    m <- roleModel(pNeut)
    
    # modify initial state so we have:
    #    - one sp close to the optim
    #    - equal abundances in metacomm 
    #    - initial sp is one not close to optim
    
    m@modelSteps[[1]]@metaComm@spAbund <- rep(1/3, 3) 
    m@modelSteps[[1]]@metaComm@spTrait <- matrix(c(-10, 0, 10))
    m@modelSteps[[1]]@localComm@indSpecies <- rep(1, 100) 
    m@modelSteps[[1]]@localComm@indTrait <- matrix(-3, nrow = 100, ncol = 1)
    
    rNeut <- runRole(m)
    
    # removing "burn-in" (1:10) from time series
    a <- ats(rNeut)[-(1:10), ]
    
    a <- apply(a, 2, function(x) {
        c(mean(x), sd(x))
    })
    
    # sp closest to optim is 2, is 2 within 95% CI of other spp?
    sp2v1 <- a[1, 2] >= a[1, 1] - 2 * a[2, 1] &
        a[1, 2] <= a[1, 1] + 2 * a[2, 1]
    sp2v3 <- a[1, 2] >= a[1, 3] - 2 * a[2, 3] &
        a[1, 2] <= a[1, 3] + 2 * a[2, 3]
    
    expect_true(sp2v1 & sp2v3)
})


test_that("species with most distinct trait is not dominant", {
    pNeut@niter <- 100000L
    pNeut@niterTimestep <- 500L
    pNeut@env_comp_delta <- 0 # as if this was comp
    m <- roleModel(pNeut)
    
    # modify initial state so we have:
    #    - one sp close to the optim
    #    - equal abundances in metacomm 
    #    - initial sp is one not close to optim
    
    m@modelSteps[[1]]@metaComm@spAbund <- rep(1/3, 3) 
    m@modelSteps[[1]]@metaComm@spTrait <- matrix(c(-10, 0, -10))
    m@modelSteps[[1]]@localComm@indSpecies <- rep(1, 100) 
    m@modelSteps[[1]]@localComm@indTrait <- matrix(-3, nrow = 100, ncol = 1)
    
    rNeut <- runRole(m)
    
    # removing "burn-in" (1:10) from time series
    a <- ats(rNeut)[-(1:10), ]
    
    a <- apply(a, 2, function(x) {
        c(mean(x), sd(x))
    })
    
    # sp closest to optim is 2, is 2 within 95% CI of other spp?
    sp2v1 <- a[1, 2] >= a[1, 1] - 2 * a[2, 1] &
        a[1, 2] <= a[1, 1] + 2 * a[2, 1]
    sp2v3 <- a[1, 2] >= a[1, 3] - 2 * a[2, 3] &
        a[1, 2] <= a[1, 3] + 2 * a[2, 3]
    
    expect_true(sp2v1 & sp2v3)
})

