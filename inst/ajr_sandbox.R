library(roleR)


neutp <- untbParams(individuals_local = 100, 
                    individuals_meta = 1000, species_meta = 100, 
                    speciation = 0.01, dispersal_prob = 0.2, 
                    init_type = 'oceanic_island', 
                    niter = 1000, niterTimestep = 10)


neut <- roleModel(neutp)
foo <- as(iterModel(neut), 'roleExperiment')
foo@experimentMeta

getSumStats(foo, list(rich = richness))





p <- roleParams(individuals_local = 100, individuals_meta = 100000, 
                species_meta = 50, speciation_local = 0.00075, speciation_meta = NULL, 
                extinction_meta = 0.05, trait_sigma = 1, env_sigma = 1, 
                comp_sigma = 1, dispersal_prob = 0.1, mutation_rate = 0.01, 
                equilib_escape = 1, alpha = 1000, num_basepairs = 500, 
                init_type = 'oceanic_island', 
                niter = 1000)




lvParams
coexistenceParams


m <- roleModel(p)

# m@modelSteps[[1]]@localComm@indSpecies <- rep(3, 100)
# m@modelSteps[[1]]@localComm@spAbund[c(1, 3)] <- c(0, 100)

foo <- iterModel(m)


ex <- as(foo, 'roleExperiment')
class(ex@modelRuns[[1]])


getSumStats(ex, list(rich = richness, hillAbund = hillAbund))



tre <- as(foo@modelSteps[[6]]@phylo, 'phylo')
plot(tre)

rawAbundance(foo@modelSteps[[2]])
foo@modelSteps[[1]]@metaComm

foo@modelSteps[[1]]@localComm@indSpecies



# ----
bla <- model@modelSteps[[1]]@localComm@indSpecies
bla


model@modelSteps[[1]]@localComm@indSpecies
