library(testthat)
library(roleR)
context('roleParams constructs and runs without error')

# works
test_that('untb constructor runs without error'){
    neutp <- untbParams(individuals_local = 1000, individuals_meta = 1000000, species_meta = 20, 
                        speciation = 0.1, dispersal_prob = 0.2, init_type = 'oceanic_island',
                        niter = 1000, niterTimestep = as.integer(100))    
}