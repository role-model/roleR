# TEST INIT
params <- list(species_meta = 100,
               individuals_meta = 10000,
               individuals_local = 1000,
               dispersal_prob = 0.1,
               speciation_local = 0.01)

#testSim
local <- localComm(abundance = numeric(), traits = matrix(numeric()), pi = numeric(), Smax = numeric())
local@abundance <- rep(0, params$species_meta * 100)
local@Smax <- params$species_meta
local@traits <- c(ape::rTraitCont(out$phylo, sigma = params$sigma_bm),
  rep(NA, params$species_meta * 99))
meta <- metaComm(abundance = numeric(), traits = matrix(numeric()), Smax = numeric())
phylo <- TreeSim::sim.bd.taxa(params$species_meta, numbsim = 1,
                              lambda = params$speciation_meta,
                              mu = params$extinction_meta)[[1]]
model <- roleModel(params, nstep = 10, nsim = 1)

# TEST BIRTH
source("R/comm.R")
source("R/birth.R")
foo <- localComm(1:10,matrix(1:100,nrow = 10, ncol = 10),1:10,10)
foo@abundance[1]
foo <- birth(foo,1)
foo@abundance[1]

source("R/rolePhylo.R")
source("R/roleSim.R")
source("R/comm.R")
p <- list(species_meta = 15,
          individuals_meta = 500,
          individuals_local = 100,
          dispersal_prob = 0.5,
          speciation_local = 0.1,
          extinction_meta = 0.8,
          speciation_meta = 1,
          trait_sigma = 0.1)

sim <- .initSim(p)
rm(sim)
sim@metaComm@traits
.iterSim(sim,10)
