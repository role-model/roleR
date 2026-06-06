# test_that("untbParams models run without error", {
#     p <- untbParams(individuals_local = 100, 
#                     individuals_meta = 1000, 
#                     species_meta = 100, 
#                     speciation = 0, 
#                     dispersal_prob = 0.1, 
#                     init_type = 'oceanic_island',
#                     niter = 300000, niterTimestep = 10000)
#     
#     expect_no_error(runRole(roleModel(p)))
# })

test_that("roleParams constructor defaults run without error", {
    p <- roleParams() 
    expect_no_error(runRole(roleModel(p)))
})

test_that("roleParams constructor returns error if timestep > niter", {
    expect_error(roleParams(niter=10, niterTimestep = 100))
})


test_that("roleParams errors out if length `timestep` or `niter` are > 1", {
    expect_error(roleParams(niter = c(10, 1), niterTimestep = c(1, 100)))
})


test_that("roleParams with a user-supplied iterfun works in runRole", {
    f <- function(i) {
        return(i / 1000)
    }
    
    p <- roleParams(speciation_local = f, niter = 100)
    
    expect_no_error(runRole(roleModel(p)))
})


test_that("roleParams with init_type = 'bridge_island' runs without error", {
    p <- roleParams(init_type = 'bridge_island', 
                    niter = 1000, niterTimestep = 10)
    m <- runRole(roleModel(p))
    
    expect_no_error(runRole(m))
})

