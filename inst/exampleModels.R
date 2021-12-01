# exampleModels contains examples of created sim models
library(roleR)
library(Rcpp)
sim <- initSim()
sim <- iterSim(sim,100)

#note - all local, meta, phylo functions take indices starting at 0 rather than 1 
sim$death()
sim$local$abundance
sim$birth()
sim$local$abundance
sim$meta$Smax
sim$immigration()
sim$meta$traits
sim$immigration()
sim$immigration()
sim$speciation()
sim$local$abundance
sim$death()
sim$birth()
sim <- iterSim(sim,1)
sim$local$abundance
