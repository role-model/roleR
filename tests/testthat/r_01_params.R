test_that("untbParams models run without error", {
    
    # bumping to 30k and 10k bugs out - verify that this is just a high speciation overflow
    p <- untbParams(individuals_local = 100, individuals_meta = 1000, 
                         species_meta = 100, 
                         speciation = 0, 
                         dispersal_prob = 0.1, init_type = 'oceanic_island',
                         niter = 30000, niterTimestep = 10000)   
    expect_error(runRole(roleModel(p)), NA)
})

test_that("roleParams constructor defaults run without error", {
    
    p <- roleParams() 
    expect_error(runRole(roleModel(p)), NA)
})

test_that("roleParams constructor returns error if timestep > niter", {
    expect_error(roleParams(niter=10,niterTimestep = 100))
})
test_that("roleParams constructor returns error if timestep or niter are >1 length or are not ints", {
    expect_error(roleParams(niter=c(10,1),niterTimestep = c(1,100)))
    expect_error(roleParams(niter=10.1,niterTimestep = 1.1))
})
    