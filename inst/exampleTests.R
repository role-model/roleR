
# clean and rebuild package first

# run this if "moving to final location" errors appear on build
Sys.setenv(R_INSTALL_STAGED = FALSE)

library("Rcpp")
library("roleR")
library("ape")

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

install.packages("devtools")
library(devtools)
install_github("ajrominger/pika")

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

#todo
#get rccp package working DONE
#finish speciation DONE WITH THINGS TO CHECK
#document objects, roleSim, & iterSimCpp DONE
#re-subclass local & meta DONE
#phylo coercion ape to role & role to c++ DONE
#add traitMax element to localComm DONE
#move function modules to objects for testing DONE
#fix bug where iterSim not updating

#move tests over
#coerce

#unable to load shared object roleR.so

#qs
# should the Cpp functions take indices starting at 1 or 0? I'm leaning towards 1, partially so they are compatible with tests and whatnot
# investigate conventions with Rcpp
# shoot an email with the package skeleton package , new repo in roleModel organization
# maximize portability & ease
# get broader group feedback about package
# testing framework k

Rcpp::Rcpp.package.skeleton()
Rcpp:package.skeleton()
