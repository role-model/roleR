# install pika
install.packages("devtools")
library(devtools)
install_github("ajrominger/pika")

# clean and rebuild package first
# Build -> Clean and Rebuild

# run this if "moving to final location" errors appear on build
Sys.setenv(R_INSTALL_STAGED = FALSE)

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

#todo
# fix speciation sampling bug
# finish testing speciation
# start working on environmental filtering
# start working on competitive filtering
#big todo
# add metadata exports
# add fitting to real data
# add intraspecific speciation 

#qs for group
# should the Cpp functions take indices starting at 1 or 0? I'm leaning towards 1, partially so they are compatible with tests and whatnot
# should localcomm, metacomm, phylocomm, etc be accessible to user? what methods should be exposed to user

