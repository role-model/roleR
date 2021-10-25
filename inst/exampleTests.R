
library("Rcpp")
library("roleR")

#------

sourceCpp("R/commCpp.cpp")
loadModule("commCpp")

m <- new(metaCommCpp, rep(1:10), matrix(), 10)
l <- new(localCommCpp, rep(1:10), matrix(), 10, rep(1:10))

#------

sourceCpp("R/paramValuesCpp.cpp")
loadModule("paramValsCpp")

pv <- new(paramValuesCpp)
pv$individuals_local

#------

sourceCpp("R/roleParamsCpp.cpp")
loadModule("paramsCpp")

p <- new(roleParamsCpp, pv, 1, "sim")
p$values$individuals_meta

#------

sourceCpp("R/rolePhyloCpp.cpp")
loadModule("phyloCpp")

phylo <- TreeSim::sim.bd.taxa(p$values$species_meta, numbsim = 1,
                            lambda = p$values$speciation_meta,
                            mu = p$values$extinction_meta, complete = FALSE)[[1]]

apeToPhyloCpp <- function(phylo){
    n = ape::Ntip(phylo)

    e = phylo$edge
    l = phylo$edge.length
    tipNames = phylo$tip.label
    tipAge = ape::node.depth.edgelength(phylo)[1:n]
    alive = rep(TRUE,n)
    alive[tipAge < max(tipAge)] = FALSE;
    scale = 1;
    out = new(rolePhyloCpp,n,e,l,alive,tipNames,scale)
    return(out)
}

phy <- apeToPhyloCpp(phylo)

#------

sourceCpp("R/roleModelCpp.cpp")
loadModule("modelCpp")

r <- new(roleModelCpp,l,m,phy,p)

sourceCpp("R/birthCpp.cpp")
loadModule("birthCpp")

l <- birthL(l,1)

r <- birthR(r)

#------

sourceCpp("src/speciationCpp.cpp")
loadModule("speciationCpp")

#------
#ROLE SIM MODEL TESTS

source("R/roleSim.R")
# initSim works
model <- .initSim()

install.packages("devtools")
library(devtools)
install_github("ajrominger/pika")

sourceCpp("R/iterSimCpp.cpp")
loadModule("iterSimCpp")

#iterSim runs, but doesnt work properly
model <- iterSim(model,100)

#todo
#get rccp package working DONE
#finish speciation DONE - TEST
#fix bug where iterSim not updating
#document objects, roleSim, & iterSimCpp DONE
#re-subclass local & meta DONE
#phylo coercion ape to role & role to c++ DONE
#add traitMax element to localComm DONE

# run this if "moving to final location" errors appear on build
Sys.setenv(R_INSTALL_STAGED = FALSE)
