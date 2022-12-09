
library(testthat)
library(roleR)
context('priors and iterfuns run without error')

# works
test_that('creating a rolePriors runs without error'){
    my_rnorm <- function(){
        return(rnorm(1),mean=0.5,sd=0.1)
    }
    params <- rolePriors(individuals_local = my_rnorm, individuals_meta = my_rnorm,
                         species_meta = my_rnorm, speciation_local = my_rnorm, #0.1
                         speciation_meta = my_rnorm, extinction_meta = my_rnorm, #1
                         trait_sigma=my_rnorm,comp_sigma = my_rnorm, env_sigma=my_rnorm,dispersal_prob = my_rnorm, mutation_rate = my_rnorm,
                         equilib_escape = my_rnorm, num_basepairs = 250,
                         init_type = 'oceanic_island', niter = 1000, niterTimestep = 100)
    # fix rep
}
