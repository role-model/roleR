 # make sure you have the current version of R and Rtools
library("roleR")

#------

m <- new(metaCommCpp, rep(1:10), matrix(), 10)

l <- new(localCommCpp, rep(1:10), matrix(), 10, rep(1:10))

#------

pv <- new(paramValuesCpp)
pv$individuals_local

#------

p <- new(roleParamsCpp, pv, "sim", 1)
p$values$individuals_meta

#------

phylo <- TreeSim::sim.bd.taxa(p$values$species_meta, numbsim = 1,
                            lambda = p$values$speciation_meta,
                            mu = p$values$extinction_meta, complete = FALSE)[[1]]

phylo <- .apeToRolePhylo(phylo)

phylo <- .rolePhyloToCpp(phylo)

#------

r <- new(roleModelCpp,l,m,phylo,p)

#------


#------
#ROLE SIM MODEL TESTS

# initSim works
model <- initSim()

model$local$birth(1)
model$local$abundance

#iterSim runs, but doesnt work properly
model <- iterSim(model,100)

library(Rcpp)
source("R/roleSim.R")
sourceCpp("src/roleModelCpp.cpp")
sourceCpp("src/modules.cpp")
sourceCpp("test_modules.cpp")

#these work
model$local$death(1)
model$local$birth(1)
model$local$immigration(1)

model$phylo$death(1)
model$phylo$speciation(1)

#these dont work as sampling isn't correct
model$birth()
model$death()
model$speciation()
model$immigration()

#benchmark
library(microbenchmark)
sim <- initSim()
microbenchmark(iterSim(sim, 100))
microbenchmark(iterSim(sim,10))

saveRDS(r, file = "model.Rd")
read <- readRDS("model.Rd")
read@localComm@abundanceS

fileConn<-file("output.txt")
writeLines(c("Hello","World"), fileConn)
close(fileConn)

#  examples of created sim model use cases (OLD)
library(roleR)
sim <- initSim()
sim2 <- sim$copy()
sim2$death()


sim$local$abundance_indv
sim2$local$abundance_indv

sim$print = FALSE
sim$local$print = FALSE
iterSim(sim,nsteps=100,50,print=FALSE)


sim$local$species_ids
sim$local$abundance_sp
View(sim$local$traitdiffs)

# traitdiffs and abundance_sp not updating properly 
sim <- iterSim(sim,1,true)
sim$local$species_ids

