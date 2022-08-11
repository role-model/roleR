reticulate::conda_install("r-reticulate",
                          packages="newick")



msprime <- reticulate::import('msprime')
newick <- reticulate::import('newick')

library(roleR)
p <- roleParams(individuals_local = 100, individuals_meta = 10000, 
                species_meta = 5, speciation_local = 0, speciation_meta = 1, 
                extinction_meta = 0.5, trait_sigma = 1, env_sigma = 1, 
                comp_sigma = 1, dispersal_prob = 0.1, mutation_rate = 0.01, 
                equilib_escape = 1, alpha = 1000, num_basepairs = 500, 
                init_type = 'oceanic_island', 
                niter = 1000, niterTimestep = 200)

mod <- roleExperiment(list(p))
mod <- iterExperiment(mod)



ex <- as(mod, 'roleExperiment')
ex@experimentMeta


foo <- getFinalState(mod)

tre <- as(foo@modelSteps[[1]]@phylo, 'phylo')
tre$edge

plot(tre)
axis(1)

treNewick <- ape::write.tree(tre)

metaSAD <- as.list(mod@params@individuals_meta * mod@modelSteps[[1]]@metaComm@spAbund)
names(metaSAD) <- tre$tip.label

metaSAD <- list(t1 = )


d <- msprime$Demography$from_species_tree(treNewick, 
                                          time_units = 'myr', 
                                          initial_size = reticulate::dict(metaSAD), 
                                          generation_time = 1, growth_rate = 0)





simGenData <- function(x) {
    # the final state of the simulation
    lastState <- getFinalState(x)
    
    # make metadata into nice table
    mdtab <- as(x, 'roleExperiment')@metadata
    
    # make a timesteps by species matrix
    x <- getSumStats(mod, list(abund = rawAbundance))
    
    l <- max(sapply(x$abund, length))
    
    timeBySpp <- lapply(x$abund, function(n) {
        c(n, rep(0, l - length(n)))
    })
    
    names(timeBySpp) <- NULL
    
    timeBySpp <- do.call(rbind, timeBySpp)
    
    # scale up to Ne
    sppNe <- timeBySpp * mdtab$alpha
    
    # take harmonic mean--that's best 
}



demography <- msprime$Demography$from_species_tree(
    "(((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)", 
    time_units = "myr",
    initial_size = 10^4, # make these harmonic means from getSumStats
    generation_time = 1, 
    growth_rate = 0)


# add pop split in a "trunk-like" topology see here:
# https://tskit.dev/msprime/docs/stable/demography.html#sec-demography-events-population-split)
# that's how you could have meta and local pops (with migration from meta added
# after the split)
# 
# if simulating multiple time steps, maybe use explicit population size setting params?
# see here: 
# https://tskit.dev/msprime/docs/stable/demography.html#population-parameters-change
# otherwise, if only simulating at end of run, then use harmonic mean...but need to record 
# that somehow...?...also only do harmonic mean from most recent 0


ts <- msprime$sim_ancestry(
    samples = list(msprime$SampleSet(6, population = 'human'), 
                   msprime$SampleSet(3, population = 'gorilla', time = 50), 
                   msprime$SampleSet(5, population = 'gorilla', time = 0)),
    demography = demography,
    ploidy = 1, 
    sequence_length = 500
)

cat(ts$draw_text())

x <- msprime$sim_mutations(ts, rate = 1e-06)

# groups of populations by time of sample
g <- split(as.integer(x$samples()), 
           paste(x$nodes_population, 
                 x$nodes_time, 
                 sep = '_')[x$nodes_flags == 1])

names(g) <- NULL

x$diversity(g)
x$as_fasta()


