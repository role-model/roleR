
test_that("model runs without error and does not end in all 1 species", {
    
    # create a test model and get data from a step and the params
    m <- quickModel()
    expect_true(length(unique(m@modelSteps[[10]]@localComm@indSpecies)) > 1)
})

test_that("big model with high speciation runs without error", {
    params <- untbParams(
        individuals_local = 900,
        individuals_meta = 9000,
        species_meta = 900,
        speciation = 0.1,
        dispersal_prob = 0.5,
        init_type = "oceanic_island",
        niter = 9000,
        niterTimestep = 10
    )
    
    m <- runRole(roleModel(params))
})

test_that("model copying in C++ does not result in all timesteps being equal", {
    
    m <- quickModel()
    #setequal(m@modelSteps[[1]]@localComm@indSpecies,m@modelSteps[[10]]@localComm@indSpecies)
    #m@modelSteps[[10]]@localComm@indSpecies
    expect_false(setequal(m@modelSteps[[1]]@localComm@indSpecies,m@modelSteps[[10]]@localComm@indSpecies))
})

test_that("testing the limits of how many iterations we can run"){
    params <- untbParams(
        individuals_local = 1000,
        individuals_meta = 10000,
        species_meta = 50,
        speciation = 0, # might break
        dispersal_prob = 0.1,
        init_type = "oceanic_island",
        niter = 100000000, 
        niterTimestep = 1000000
    )
    
    m <- runRole(roleModel(params))
}

test_that("comparing runtime of R implementation in a simple example"){
    x <- proc.time()
    params <- untbParams(
        individuals_local = 1000,
        individuals_meta = 10000,
        species_meta = 50,
        speciation = 0, # might break
        dispersal_prob = 0.1,
        init_type = "oceanic_island",
        niter = 1000000, 
        niterTimestep = 100000
    )
    m <- runRole(roleModel(params))
    proc.time() - x
    # takes 17.5 seconds relative to 45 in base R with missing functionality
}

test_that("when a model is run the supplied model is NOT modified in place"){
    m <- quickModelNonRun()
    mrun <- runRole(m)
    
    expect_true(is.null(m@modelSteps[[2]]))
    expect_false(is.null(mrun@modelSteps[[2]]))
}