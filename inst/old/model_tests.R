
# fully ran 1 time
library(roleR)
p <- FALSE
r <- dummyModel(R=FALSE,run=FALSE)

for(i in 1:1000)
{
  print(i)
  #iterate 1 time
  r <- iterSim(r,1,1,FALSE)
  d <- new(roleDataCpp,r$local,r$meta,r$phylo)
  d <- roleDataFromCpp(d)
  saveRDS(d,"test_data.roledata")
  Sys.sleep(0.4)
}
# for a while was going exactly to 497-499 
library(pryr)
mem_used()
object.size(r)
test <- readRDS("test_data.roledata")
r$local$abundance_indv <- test@localComm@abundanceIndv
r$local$species_ids <- test@localComm@speciesIDsIndv
r$local$traits <- test@localComm@traitsIndv
r$local$abundance_sp <- test@localComm@abundanceSp
r$local$traits_sp <- test@localComm@traitsSp
r$meta$abundance <- test@metaComm@abundanceSp
r$meta$traits <- test@metaComm@traitsSp
r$phylo$n <- test@phylo@ntips
r$phylo$e <- test@phylo@edges
r$phylo$l <- test@phylo@lengths
r$phylo$alive <- test@phylo@alive
r$phylo$tipNames <- test@phylo@tipNames
r$phylo$scale <- test@phylo@scale
  
r$params$dispersal_prob
validObject(p)
print(Sys.time() - start_time)

params <- roleParams(nruns=10,niter=1000,defaults = TRUE)
sim <- roleSim(params)

start_time = Sys.time()
r <- dummyModel(R=FALSE)
iterSim(r,1000,100,p)
print(Sys.time() - start_time)

iterSim(r,1000,100,p)
r$local$species_ids[114]
r$local$traits
r$local$traits
r$local$traits_sp
r$local$traits[r$local$species_ids == 14]
r$local$species_ids
library(roleR)
r <- scratchModel()
r$print = FALSE
iterSim(r,1000,100,TRUE)

test <- readRDS("test.roledata")
test@metaComm
r <- scratchModel()
r$local$traits_sp

# traits sp is getting super big or negative
test@localComm@traitsSp

test@localComm@abundanceSp
# causes traits indv to get negative too 

#abundance sp becomes proportional 

test

# roleParams object
params <- roleParams(nrun=1,niter=100,niter_timestep=10,defaults=TRUE)
params@values[[1]]
params <- setDefaultParams(params) # set defaults works
params <- setParam(params,"species_meta",15) # set param works
cparams <- stretchAndSampleParams(params)
cparams@values[[1]]

# initSim seemingly works
sim <- initSim(cparams@values[[1]],type="bridge_island")

# roleSim
sim <- roleSim(params,print=TRUE)
out <- iterSim(sim,1,100,10)

cparams <- stretchAndSampleParams(params)
parlist <- cparams@values[[1]]
parlist[["species_meta"]]
parlist[["comp_sigma"]]

length(params@priors)
class(sim)
iterSim(sim,100,10, TRUE)

sim$
iterSim(sim, 1,100,10)

#all
library(roleR)
i = 10
params <- roleParams(nrun=1,niter=i,niter_timestep=10,defaults=TRUE)
cparams <- stretchAndSampleParams(params)
parlist <- cparams@values[[1]]
parlist
sim <- initSim(parlist,type="bridge_island")
sim$phylo$n
sim$local$Smax
sim$meta$Smax
log <- TRUE
sim$print <- log
sim$local$print <- log
iterSim(sim,1,10,log)
object.size(sim) 

sim$phylo$n
sim$local$Smax
sim$meta$Smax

sim <- roleModelFromCpp(sim)

#all2
library(roleR)
params <- roleParams(nrun=1,niter=1000,niterTimestep=10,defaults=TRUE)
sim <- roleSim(params,print=FALSE)
cparams <- stretchAndSampleParams(params)
parlist <- cparams@values[[1]]
sim <- initSim(parlist,type="bridge_island")
iterSim(sim,100,10,FALSE)

#all3
library(roleR)
i <- 10000
params <- roleParams(nrun=1,niter=i,niterTimestep=50,defaults=TRUE)
params <- setParam(params,"speciation_local",0.75)
params <- setParam(params,"speciation_meta",0.75)
cparams <- stretchAndSampleParams(params)
parlist <- cparams@values[[1]]
sim <- initSim(parlist,type="bridge_island",niter=i)

iterSim(sim,i,5,FALSE)

sim$local$abundance_sp
sim$phylo$n
sim$local$Smax
sim$meta$Smax

sim <- roleModelFromCpp(sim)

10 %% 100
# loop test
library(roleR)
r <- scratchModel()
d <- new(roleDataCpp,r$local,r$meta,r$phylo)
r$timeseries[[1]] <- d
rr <- roleModelFromCpp(r)

d <- new(roleDataCpp,r$local,r$meta,r$phylo)
r$phylo$
r$local$traits_sp
rd <- roleDataFromCpp(d)
rd@
rr <- roleModelFromCpp(r)
length(r$timeseries)
d <- roleDataFromCpp(r)

library(roleR)
r <- scratchModel()
r$timeseries[[1]] <- r$copyData(0)
r$print = TRUE
for(i in 1:1000)
{
  print(i)
  #iterate 1 time
  iterSim(r,1,100,TRUE)
  d <- new(roleDataCpp,r$local,r$meta,r$phylo)
  #rr <- roleModelFromCpp(r)
  #saveRDS(rr,"test.rolemodel")
  #roleDataFromCpp(r$timeseries[0])
}

