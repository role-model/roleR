test_that("experiment runs two models with different params without error", {
    
    pu <- untbParams(individuals_local = 100, individuals_meta = 1000, 
                    species_meta = 50, 
                    speciation = 0.2, 
                    dispersal_prob = 0.1, init_type = 'oceanic_island',
                    niter = 1000, niterTimestep = 100)   
    p <- roleParams(individuals_local = 100, individuals_meta = 1000,
                    species_meta = 50, speciation_local = 0.05, 
                    speciation_meta = 0.05, extinction_meta = 0.05, env_sigma = 0.5,
                    neut_delta = 0, env_comp_delta=0.5, alpha=1,
                    trait_sigma=1,comp_sigma = 0.5, dispersal_prob = 0.1, mutation_rate = 0.01,
                    equilib_escape = 1, num_basepairs = 250,
                    init_type = 'oceanic_island', niter = 1000, niterTimestep = 100)
    
    exp <- runRole(roleExperiment(list(p,pu)))
})