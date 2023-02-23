
# function to create a model using a set of reasonable params
# used in testthat tests
quickModel <- function(){
    p <- roleParams(individuals_local = 100, individuals_meta = 1000,
                    species_meta = 10, speciation_local = 0.1, 
                    speciation_meta = 0.1, extinction_meta = 0.05, env_sigma = 0.5,
                    trait_sigma=1,comp_sigma = 0.5, dispersal_prob = 0.1, mutation_rate = 0.01,
                    equilib_escape = 1, num_basepairs = 250,
                    init_type = 'oceanic_island', niter = 100, niterTimestep = 10)
    
    return(runRoLE(roleModel(p))) 
}