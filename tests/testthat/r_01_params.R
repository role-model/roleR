test_that("untbParams models run without error", {
    
    p <- untbParams(individuals_local = 100, individuals_meta = 1000, 
                         species_meta = 50, 
                         speciation = 0.2, 
                         dispersal_prob = 0.1, init_type = 'oceanic_island',
                         niter = 3000, niterTimestep = 1000)   
    expect_error(runRoLE(roleModel(p)), NA)
})