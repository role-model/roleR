
library(roleR)
params <- roleParams(individuals_local = 100, individuals_meta = 1000,
                     species_meta = 500, speciation_local = 0.1, #0.1
                     speciation_meta = 1, extinction_meta = 0.8, #1
                     trait_sigma=1,comp_sigma = 0.1, env_sigma=0.1,dispersal_prob = 0.1, mutation_rate = 0.01,
                     equilib_escape = 1, num_basepairs = 250,
                     init_type = 'oceanic_island', niter = 1000, niterTimestep = 1000)

params <- roleParams(niter=100)
params@individuals_local




model <- roleModel(params)

data_list <- iterModel(model,F)
state <- getFinalState(data_list)
state@
state@
foo <- as(data_list, 'roleExperiment')
boo <- getSumStats(foo, list(abund = rawAbundance, rich = richness))


# note - probably have to add one to indSpecies afterward
data_list[[1]]@localComm@indSpecies
data_list[[10]]@localComm@indSppTrt
prepRebuild()

p2 <- params
p3 <- params

install.packages("parallel")
params_list <- list(params,p2,p3)
library(parallel)

experiment <- roleExperiment(params_list, priors_list) # throws "attempt to replicate an object of type S4" error
runs <- experiment@modelRuns
experiment <- iterExperiment(exp)


install.packages("snow")
z=vector('list',4)
z=1:4
system.time(lapply(z,function(x) Sys.sleep(1)))
cl<-makeCluster(2,type="SOCK")
system.time(clusterApply(cl, z,function(x) Sys.sleep(1)))
stopCluster(cl)

clusterApply(cl,)
