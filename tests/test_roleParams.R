library(testthat)
library(roleR)
context('roleParams constructs and runs without error')

# works
test_that('minimum constructor runs without error'){
    params <- roleParams(individuals_local = 100, individuals_meta = 500,
                         species_meta = 500, speciation_local = 0.1, 
                         speciation_meta = 1, extinction_meta = 0.8, #env_sigma = 0.5,
                         trait_sigma=1,comp_sigma = 0.1,# dispersal_prob = 0.1
                         init_type = 'oceanic_island', niter = 1000, niterTimestep = 10)
}

# works
test_that('untb constructor runs without error'){
    neutp <- untbParams(individuals_local = 1000, individuals_meta = 1000000, species_meta = 20, 
                        speciation = 0.1, dispersal_prob = 0.2, init_type = 'oceanic_island',
                        niter = 1000, niterTimestep = as.integer(100))    
}

# works
# NOTE - if speciation_meta is less than 0.7? putting it into a model hangs
# NOTE - if speciation_local is high an index is exceeded Index out of bounds: [index=5050; extent=5050]
# probably something in the phylo
# Index out of bounds: [index=10000; extent=10000].
test_that('untb constructor runs without error'){
    
    params <- untbParams(individuals_local = 100, individuals_meta = 1000, 
                        species_meta = 50, 
                        speciation = 0.2, 
                        dispersal_prob = 0.1, init_type = 'oceanic_island',
                        niter = 30000, niterTimestep = 1000)   
    
    model <- roleModel(neutp)
    model <- runRoLE(model)
}