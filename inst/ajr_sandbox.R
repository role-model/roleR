library(roleR)

p <- roleParams(individuals_local = 100, individuals_meta = 10000, 
                species_meta = 5, speciation_local = 0.005, speciation_meta = 1, 
                extinction_meta = 0.5, trait_sigma = 1, env_sigma = 1, 
                comp_sigma = 1, dispersal_prob = 0.1, mutation_rate = 0.01, 
                equilib_escape = 1, alpha = 1000, num_basepairs = 500, 
                init_type = 'oceanic_island', 
                niter = 1000, niterTimestep = 200)

m <- roleModel(p)

foo <- iterModel(m)

foo@modelSteps[[2]]@phylo@e[1:11, ]

tre <- as(foo@modelSteps[[2]]@phylo, 'phylo')
plot(tre)



foo@modelSteps[[1]]@localComm@indSpecies



# ----
bla <- model@modelSteps[[1]]@localComm@indSpecies
bla


model@modelSteps[[1]]@localComm@indSpecies
