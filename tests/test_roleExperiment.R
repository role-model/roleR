library(testthat)
library(roleR)
context('roleExperiment constructs and runs without error')

# works
test_that('experiment constructor runs without error'){
    params <- roleParams(individuals_local = 100, individuals_meta = 1000,
                         species_meta = 500, speciation_local = 0.1, 
                         speciation_meta = 1, extinction_meta = 0.8, 
                         trait_sigma=1,comp_sigma = 1, dispersal_prob = 0.1, mutation_rate = 0.01,
                         equilib_escape = 1, num_basepairs = 250,
                         init_type = 'oceanic_island', niter = 10, niterTimestep = 10) 
    p2 <- params
    p3 <- params
    params_list <- list(params,p2,p3)
    experiment <- roleExperiment(params_list)
}

# works
test_that('experiment models iterated over successfully, unparallelized'){
    params <- roleParams(individuals_local = 100, individuals_meta = 1000,
                         species_meta = 500, speciation_local = 0.1, 
                         speciation_meta = 1, extinction_meta = 0.8, 
                         trait_sigma=1,comp_sigma = 1, dispersal_prob = 0.1, mutation_rate = 0.01,
                         equilib_escape = 1, num_basepairs = 250,
                         init_type = 'oceanic_island', niter = 10, niterTimestep = 10) 
    p2 <- params
    p3 <- params
    params_list <- list(params,p2,p3)
    experiment <- roleExperiment(params_list)
    experiment@modelRuns[[1]]@modelSteps[[2]]
    experiment <- iterExperiment(experiment)
    experiment@modelRuns[[1]]@modelSteps[[2]]
}

# works
test_that('experiment models iterated over successfully, parallelized'){
    params <- roleParams(individuals_local = 100, individuals_meta = 1000,
                         species_meta = 500, speciation_local = 0.1, 
                         speciation_meta = 1, extinction_meta = 0.8, 
                         trait_sigma=1,comp_sigma = 1, dispersal_prob = 0.1, mutation_rate = 0.01,
                         equilib_escape = 1, num_basepairs = 250,
                         init_type = 'oceanic_island', niter = 10000000, niterTimestep = 100000) 
    p2 <- params
    p3 <- params
    params_list <- list(params,p2,p3)
    experiment <- roleExperiment(params_list)
    experiment@modelRuns[[1]]@modelSteps[[100]]
    experiment <- iterExperiment(experiment,cores=3) #1 = 26 seconds, 3 = 16 seconds
    experiment@modelRuns[[1]]@modelSteps[[100]]
}