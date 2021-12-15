# exampleModels contains examples of created sim models
library(roleR)
library(Rcpp)

sim <- initSim()

sim <- iterSim(sim,1000)

meme <- matrix(0, nrow = 10000,ncol= 10000)
View(meme)
