
# function to create a model using a set of reasonable params
# used in testthat tests
# will eventually no longer be exported
quickModel <- function(){
    # its for some reason using the value of the last supplied non function in every fun 
    p <- roleParams(individuals_local = 100, individuals_meta = 1000, species_meta = 100, 
                     speciation_local = 0.1, speciation_meta = 0.1, extinction_meta = 0.05, dispersal_prob = 0.1,
                     
                     trait_sigma=1, env_sigma=1, comp_sigma = 0.5, neut_delta=1, env_comp_delta=1,
                     mutation_rate=0,equilib_escape = 1, alpha=50, num_basepairs = 250,
                     init_type = 'oceanic_island', niter = 100, niterTimestep = 10)
    m <- roleModel(p)
    
    return(runRole(m)) 
}

# function to create a model using a set of reasonable params
# used in testthat tests
# will eventually no longer be exported
quickExp <- function(){
    # its for some reason using the value of the last supplied non function in every fun 
    p <- roleParams(individuals_local = 100, individuals_meta = 1000, species_meta = 100, 
                    speciation_local = 0.1, speciation_meta = 0.1, extinction_meta = 0.05, dispersal_prob = 0.1,
                    trait_sigma=1, env_sigma=1, comp_sigma = 0.5, neut_delta=1, env_comp_delta=1,
                    mutation_rate=0,equilib_escape = 1, alpha=50, num_basepairs = 250,
                    init_type = 'oceanic_island', niter = 100, niterTimestep = 10)
    m <- roleModel(p)
    
    return(runRole(m)) 
}