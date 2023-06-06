# default params
p <- roleParams()

# initialized model with default params
m <- try(roleModel(p), silent = TRUE)


# tests of initialized model ----
test_that("model using default params initializes without error", {
    expect_true(inherits(m, 'roleModel'))
})


test_that("model info data.frame is of correct dimension", {
    df <- m@info
    expect_equal(nrow(df), p@niter / p@niterTimestep + 1)
})

# keep this for later
initInfo <- m@info

test_that("`modelSteps` output of model is initialized correctly", {
    msteps <- m@modelSteps
    expect_equal(length(msteps), nrow(m@info))
    expect_true(inherits(msteps[[1]], 'roleData'))
    
    notInit <- sapply(msteps[-1], is.null)
    expect_true(all(notInit))
})


# should continue to add tests to make sure sub-objects of data are valid

# tests on run model ----


# run model
m <- try(runRole(m), silent = TRUE)


test_that("model runs without error", {
    expect_true(inherits(m, 'roleModel'))
})

test_that("model info is not altered by running model", {
    expect_equal(m@info, initInfo)
})

test_that("`modelSteps` output is correctly updated by run", {
    msteps <- m@modelSteps
    expect_equal(length(msteps), nrow(m@info))
    
    allSteps <- sapply(msteps, inherits, what = 'roleData')
    expect_true(all(allSteps))
})

test_that("when a model is run the supplied model is NOT modified in place", {
    m <- roleModel(p)
    mrun <- runRole(m)
    
    expect_true(is.null(m@modelSteps[[2]]) & !is.null(mrun@modelSteps[[2]]))
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
    
    expect_true(length(m@modelSteps) == 901)
})

test_that("model copying in C++ does not result in all timesteps being equal", {
    p <- roleParams(niter = 200, speciation_local = .5) 
    # create a test model and get data from a step and the params
    m <- roleModel(p)
    m <- runRole(m)    
    
    expect_false(setequal(m@modelSteps[[2]]@localComm@indSpecies,
                          m@modelSteps[[11]]@localComm@indSpecies))
})



test_that("metacommunity logseries initialization is correct", {
    S <- 20
    N <- 1000
    s <- meteR::meteDist2Rank(meteR::sad(meteR::meteESF(S0 = S, N0 = N)))
    s <- s / sum(s)
    m <- roleR:::.lseriesFromSN(S, N)
    expect_lt(sum((s - m)^2), 0.0001)
})
