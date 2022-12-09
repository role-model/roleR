
library(testthat)
library(roleR)
context('get sum stats works')

test_that('get sum stats runs without error'){
    params <- roleParams(individuals_local = 100, individuals_meta = 1000,
                         species_meta = 500, speciation_local = 0.1, #0.1
                         speciation_meta = 1, extinction_meta = 0.8, #1
                         trait_sigma=1,comp_sigma = 0.1, env_sigma=0.1,dispersal_prob = 0.1, mutation_rate = 0.01, # comp sigma =1 causes probs to be bad
                         equilib_escape = 1, num_basepairs = 250,
                         init_type = 'oceanic_island', niter = 1000, niterTimestep = 1000)
    model <- roleModel(params)
    model <- iterModel(model,F)
    
    #still need to ask what this should do exactly
}