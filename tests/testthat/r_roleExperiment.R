# three near-default params
pp <- list(roleParams(individuals_local = 50, niter = 8), 
           roleParams(individuals_local = 75),
           roleParams(individuals_local = 100))

# initialized experiment
e <- try(roleExperiment(pp), silent = TRUE)


# tests of initialized model ----
test_that("experiment using near default params initializes without error", {
    expect_true(inherits(e, 'roleExperiment'))
})

# initialized `info` df
initInfo <- e@info

# num timesteps in each model
allN <- sapply(pp, function(p) {
    p@niter / p@niterTimestep + 1
})

test_that("experiment info data.frame is of correct dimension", {
    expect_equal(nrow(initInfo), sum(allN))
})


test_that("`modelRuns` output of model is initialized correctly", {
    # model runs is correct length
    msteps <- e@modelRuns
    expect_equal(length(msteps), nrow(e@info))
    
    # time step 0 for each model has initialized roleData
    ii <- which(e@info$timestep == 0)
    x0 <- sapply(msteps[ii], inherits, what = 'roleData')
    expect_true(all(x0))
    
    # non-initialized 
    notInit <- sapply(msteps[-ii], is.null)
    expect_true(all(notInit))
})

# include a test for the `inits` slot

# should continue to add tests to make sure sub-objects of data are valid


# tests on run experiment ----

# run experiment
erun <- try(runRole(e), silent = TRUE)


test_that("experiment runs without error", {
    expect_true(inherits(erun, 'roleExperiment'))
})

test_that("model info is not altered by running model", {
    expect_equal(erun@info, initInfo)
})

test_that("`modelRuns` output is correctly updated by run", {
    msteps <- erun@modelRuns
    expect_equal(length(msteps), nrow(erun@info))
    
    allSteps <- sapply(msteps, inherits, what = 'roleData')
    expect_true(all(allSteps))
})

test_that("when an experiment is run the supplied experiment is NOT modified in place", {
    expect_true(is.null(e@modelRuns[[2]]) & !is.null(erun@modelRuns[[2]]))
})
