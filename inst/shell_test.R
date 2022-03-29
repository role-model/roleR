# script to run model on command line for testing

library(roleR)
params <- roleParams(nrun=1,niter=100,niterTimestep=100,defaults=TRUE)
cparams <- stretchAndSampleParams(params)
parlist <- cparams@values[[1]]
model <- initSim(parlist,type="bridge_island",niter=1000)
model$print <- FALSE
model$local$print <- FALSE

print(model$local$abundance_sp)

iterSim(model,100,params@niterTimestep,FALSE)

print(model$local$abundance_sp)

