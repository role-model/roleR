setwd("roleRcpp")

library("Rcpp")

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

#ROLE SIM TESTS

source("R/roleSim.R")
.initSim()
