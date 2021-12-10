# exampleModels contains examples of created sim models
library(roleR)
library(Rcpp)

sim <- initSim()
sim <- iterSim(sim,100)
