
library(testthat)
library(roleR)
context('roleModel constructs and runs without error')

# works
test_that('model constructor runs without error'){
    params <- roleParams(individuals_local = 100, individuals_meta = 1000,
                         species_meta = 500, speciation_local = 0.1, 
                         speciation_meta = 1, extinction_meta = 0.8, 
                         trait_sigma=1,comp_sigma = 1, dispersal_prob = 0.1, mutation_rate = 0.01,
                         equilib_escape = 1, num_basepairs = 250,
                         init_type = 'oceanic_island', niter = 10, niterTimestep = 10) 
    model <- roleModel(params)
}

# works
test_that('model of 100 iters runs without crash'){
    params <- roleParams(individuals_local = 100, individuals_meta = 1000,
                         species_meta = 500, speciation_local = 0.1, 
                         speciation_meta = 1, extinction_meta = 0.8, 
                         trait_sigma=1,comp_sigma = 1, dispersal_prob = 0.1, mutation_rate = 0.01,
                         equilib_escape = 1, num_basepairs = 250,
                         init_type = 'oceanic_island', niter = 100, niterTimestep = 10) 
    model <- roleModel(params)
    model <- iterModel(params)
}
