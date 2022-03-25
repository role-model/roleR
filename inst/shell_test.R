# test script to run model on R cmd

library(roleR)
params <- roleParams(nrun=1,niter=1000,niterTimestep=100,defaults=TRUE)
cparams <- stretchAndSampleParams(params)
parlist <- cparams@values[[1]]
model <- initSim(parlist,type="bridge_island",niter=1000)

model$print <- FALSE
model$local$print <- FALSE

iterSim(model,100,params@niterTimestep,FALSE)
model$local$abundance_sp

d <- roleModelFromCpp(model)
model$timeseries[1]
saveRDS(d,"test.roledata")
model$local$abundance_sp


for(i in 1:10)
{
  print(i)
  #iterate 1 time
  iterSim(model,100,100,FALSE)
  #d <- new(roleDataCpp,r$local,r$meta,r$phylo)
  #d <- roleDataFromCpp(d)
  #saveRDS(d,"test.roledata")
}

model$