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
