library(roleR)
p <- roleParams(individuals_local = 100, individuals_meta = 1000,
                species_meta = 10, speciation_local = 0.1, 
                speciation_meta = 0.1, extinction_meta = 0.05, env_sigma = 0.5,
                trait_sigma=1,comp_sigma = 0.5, dispersal_prob = 0.1, mutation_rate = 0.01,
                equilib_escape = 1, num_basepairs = 250,
                init_type = 'oceanic_island', niter = 100, niterTimestep = 10)

model <- runRoLE(roleModel(p))

model@modelSteps[[11]]@phylo@alive
phy <- as(model@modelSteps[[11]]@phylo, 'phylo')
#phy <- ape::rphylo(10, p@speciation_meta, p@extinction_meta,T0=0)

plot(phy,use.edge.length = TRUE)
add.scale.bar()
axisPhylo()
edgelabels(phy$edge.length,cex=0.5)

#phy$edge.length <- rep(0,length(phy$edge.length))

nwk <- ape::write.tree(phy)

d <- msprime$Demography$from_species_tree(nwk, 
                                          time_units = 'gen', 
                                          initial_size = 1100, growth_rate = 0)

for 
#must find all extinct populations
d
#crap
phy$edge.length <- round(phy$edge.length)
phy$edge.length
tip.height(phy)
node.height(phy)
node.depth.edgelength(phy)
phy$tip.label
node.depth(phy)

model@modelSteps[[10]]@phylo@alive[1:30]

