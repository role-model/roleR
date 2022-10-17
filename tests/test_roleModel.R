
library(testthat)
library(roleR)
context('roleModel constructs and runs without error')

# works
test_that('model constructor runs without error'){
    params <- roleParams(individuals_local = 100, individuals_meta = 1000,
                         species_meta = 50, speciation_local = 0.3, 
                         speciation_meta = 1, extinction_meta = 0.8, env_sigma = 0.5,
                         trait_sigma=1,comp_sigma = 0.1, dispersal_prob = 0.1, mutation_rate = 0.01,
                         equilib_escape = 1, num_basepairs = 250,
                         init_type = 'oceanic_island', niter = 3, niterTimestep = 3) 
    
    model <- roleModel(p)
    model@modelSteps[[1]]@phylo@e
    model@modelSteps[[1]]@localComm@indSpecies
    foo <- roleR::runRoLE(model)
    getSumStats(foo, list(rich = richness))
    
    foo@modelSteps[[1]]@phylo@e
    foo@modelSteps[[2]]@localComm@indSpecies
    model@modelSteps[[1]]@localComm@indSpecies
    prepRebuild()
}

# works
test_that('model of 100 iters runs without crash'){
    params <- roleParams(individuals_local = 100, individuals_meta = 500,
                         species_meta = 500, speciation_local = 0.1, 
                         speciation_meta = 1, extinction_meta = 0.8, env_sigma = 0.5,
                         trait_sigma=1,comp_sigma = 0.1, dispersal_prob = 0.1, mutation_rate = 0.01,
                         equilib_escape = 1, num_basepairs = 250,
                         init_type = 'oceanic_island', niter = 1000, niterTimestep = 10)
    model <- roleModel(params)
    model <- iterModel(model,F)
}

test_that('model outputs are reasonable using untb'){
    
    neutp <- untbParams(individuals_local = 100, 
                        individuals_meta = 10000, species_meta = 50, 
                        speciation = 0.01, dispersal_prob = 0.1, 
                        init_type = 'oceanic_island', 
                        niter = 1000, niterTimestep = 10)
    neut <- roleModel(neutp)
    neut <- iterModel(neut)
    exp <- as(neut, 'roleExperiment')
    
    getSumStats(exp, list(rich = richness))
}

# NOTE - if neut delta is not specified it defaults to NA and causes "probs must be finite" error 
# NOTE - if you iterate a model twice everything gets completely screwed up 
test_that('model outputs are reasonable using normal params'){
    params <- roleParams(individuals_local = 100, individuals_meta = 10000,
                         species_meta = 50, speciation_local = 0.01, 
                         speciation_meta = 0.8, extinction_meta = 0.05, env_sigma = 1,
                         trait_sigma=1,comp_sigma = 1, dispersal_prob = 0.1, mutation_rate = 0.01,
                         equilib_escape = 1, num_basepairs = 250,
                         init_type = 'oceanic_island', niter = 1000, niterTimestep = 10, neut_delta=0,env_comp_delta = 0.8)
    model <- roleModel(params)
    #model@params@neut_delta
    model <- iterModel(model)
    exp <- as(model, 'roleExperiment')
    
    getSumStats(exp, list(rich = richness))
}
