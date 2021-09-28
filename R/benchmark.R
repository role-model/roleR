library(microbenchmark)
#hacky way of doing this for now

#default set of plausible params
p <- list(species_meta = 15,
          individuals_meta = 500,
          individuals_local = 100,
          dispersal_prob = 0.5,
          speciation_local = 0.1,
          extinction_meta = 0.8,
          speciation_meta = 1,
          trait_sigma = 0.1)
p <- roleParams(p,"sim")

#bench phylo simulation
microbenchmark(phy <- TreeSim::sim.bd.taxa(p@params$species_meta, numbsim = 1,
                            lambda = p@params$speciation_meta,
                            mu = p@params$extinction_meta, complete = FALSE)[[1]])


# bench metacommunity SAD simulation
microbenchmark(meta@abundance <- .lseriesFromSN(p@params$species_meta,
                                 p@params$individuals_meta))

#prep meta
meta <- metaComm(abundance = numeric(), traits = matrix(numeric()),
                 Smax = 0)
meta@Smax <- p@params$species_meta

# bench metacomm trait simulation
microbenchmark(meta@traits <- cbind(1:meta@Smax,
                     ape::rTraitCont(phy, sigma = p@params$trait_sigma)))

# bench iterSim
microbenchmark(.iterSim(defaultRoleModel, 10))
