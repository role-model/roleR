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
test_that("roleParams with a user-supplied iterfun works in runRole", {
    
    f <- function(i) {
        return(i / 1000)
    }
    
    p <- roleParams(speciation_local = f)
    expect_error(runRole(roleModel(p)),NA)
})

test_that("roleParams with init_type = 'bridge_island' runs a model without error", {
    
    p <- roleParams(init_type = 'bridge_island', niter = 1000, niterTimestep = 10)
    m <- runRole(roleModel(p))
    expect_error(runRole(roleModel(p)),NA)
})

test_that("roleParams with very high speciation rate runs a model without error", {
    
    # bug scales with niter and speciation rate
    p <- roleParams(speciation_local=0.99, niter=1000)
    expect_error(runRole(roleModel(p)),NA)
})

test_that("roleParams with neut delta=1 runs", {
        p <- roleParams(neut_delta=1)
    expect_error(runRole(roleModel(p)),NA)
})

# fails!!
test_that("roleParams with neut delta >0 <1  runs", {
        p <- roleParams(neut_delta=0.5)
    expect_error(runRole(roleModel(p)),NA)
})
