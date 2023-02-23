
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
    
    m <- runRoLE(roleModel(params))
})

test_that("model copying in C++ does not result in all timesteps being equal", {
    
    m <- quickModel()
    #setequal(m@modelSteps[[1]]@localComm@indSpecies,m@modelSteps[[10]]@localComm@indSpecies)
    #m@modelSteps[[10]]@localComm@indSpecies
    expect_false(setequal(m@modelSteps[[1]]@localComm@indSpecies,m@modelSteps[[10]]@localComm@indSpecies))
})