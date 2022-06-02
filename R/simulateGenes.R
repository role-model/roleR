# install miniconda, r-reticulate, and msprime
# execute on package load for next week
installGeneSimReqs <- function()
{
  library(reticulate)
  install_miniconda()
  conda_install("r-reticulate",
                packages="msprime",
                channel="conda-forge")
}

# <- dummyrole()
#run <- dummyModel(R=TRUE)
#params <- roleParams(nrun=1,niter=1000,niterTimestep=10,defaults=TRUE)
#params@values[1][[1]]$mutation_rate <- 0.1
#params@values[1][[1]]$num_basepairs <- 500
#role <- readRDS("data/example_out_in/test.roleexperiment")
# role@params@values[[1]]$num_basepairs <- 10
# error Error in py_call_impl(callable, dots$args, dots$keywords) : 
# msprime._msprime.LibraryError: The simulation model supplied resulted in a parent node having a time value <= to its child. This can occur either as a result of multiple bottlenecks happening at the same time or because of numerical imprecision with very small population sizes.
# simulateSpeciesGeneticDiv(role)

simulateSpeciesGeneticDiv <- function(role){
  
  # import msprime into R 
  library(reticulate)
  msprime <- import("msprime")
  
  # get params from role experiment
  params <- role@params
  niter <- params@niter
  niter_timestep <- params@niterTimestep
  indv_local <- params@values[[1]][["individuals_local"]]
  alpha <- params@values[["alpha"]]
  
  # for every model run in the role experiment
  for(i in 1:length(role@modelRuns))
  {
    # i = 1
    # get model and num_basepairs
    model <- role@modelRuns[[i]]
    num_basepairs <- params@values[[i]]$num_basepairs
    
    # get the model timeseries of species abundances
    abundance_ts <- getTimeseries(model,type="model_value",valueName="abundanceSp")
    # convert to a matrix 
    abundance_ts <- t(matrix(unlist(abundance_ts),ncol=niter_timestep))
    
    # delete cols where every value is 0 - don't understand why this is again, ask 
    # resolves problem of resolution where a species could exist then go extinct within one interval
    # abundance_ts <- abundance_ts[, colSums(abundance_ts) > 0]
    
    sp_indices <- which(colSums(abundance_ts) > 0)
  
    # initialize g_diversities
    g_diversities <- numeric(ncol(abundance_ts))
    # place NAs for species where no individuals are recorded to have existed 
    g_diversities[-sp_indices] <- NA
    
    # for each species (each column in the abundance timeseries matrix) 
    for(sp_index in sp_indices)
    {
      # return null to columns with only zeroes 
      # sp_index = 1
      # at a minimum save every x generations 
      # get vector of abundances over time for that species
      abundance_vect <- abundance_ts[,sp_index] * diem_scalar
      
      # get index of first nonzero abundance - this is the timestep at which the species started
      start_timestep <- min(which(abundance_vect != 0))
      end_timestep <- max(which(abundance_vect != 0))
      
      # start - end * iters per timestep
      # 1 generation = J/2 timesteps 
      # J
      # diem size - scalar for different between observed and actual counts 
      n_steps_lived <- ((end_timestep - start_timestep) * niter_timestep) / (J/2)
      
      local_pop_start <- abundance_vect[start_timestep]
      local_pop_end <- abundance_vect[end_timestep]
      
      # if the species originated in local and thus did not come from meta
      if(sp_index > length(model@timeseries[[1]]@metaComm@abundanceSp))
      {
        # sample size is species abundance at present
        # Ne is effective populatin size - harmonic mean of all abundances of all timesteps
        # Ne is alpha times harmonic mean 
        # alpha gets estimated pretty accurately
        # below no meta
        
        # num_basepairs = 10
        ts <- msprime$simulate(sample_size=as.integer(10),
                               Ne=as.integer(1e5),
                               length=as.integer(num_basepairs),
                               mutation_rate=as.integer(1e-7))
        
        bp_diversity <- ts$diversity()/num_basepairs
      }
      
      # otherwise the species immigrated from meta
      else
      {
        meta_pop_startend <- model@timeseries[[1]]@metaComm@abundanceSp[sp_index] * 100 # * individuals meta
        
        # create the meta community
        pop_meta <- msprime$PopulationConfiguration(sample_size=as.integer(meta_pop_startend),
                                                    initial_size=as.integer(meta_pop_startend))
        pop_local <- msprime$PopulationConfiguration(sample_size=as.integer(local_pop_end),
                                                     initial_size=as.integer(local_pop_start))
        
        # specify migration from local to metacommunity (backwards in time)
        mig_matrix <- matrix(c(0,params@values[1][[1]]$dispersal_prob,0,0),nrow=2,ncol=2)
        
        # disperal_rate may be a scalar of the prob
        # same for iters
        # tdiv <- niters * 1/J (J is # local individuals of species at niter)
        # moran process to wright fisher process with overlapping gens 
        
        # divergence time is equal to immigration time
        tdiv <- start_timestep
        split_event <- msprime$MassMigration(time=as.integer(tdiv),
                                             source=as.integer(1),
                                             destination=as.integer(0),
                                             proportion=as.integer(1))
        
        # params@values[1][[1]]$mutation_rate <- 0.1
        # Length of sequence to simulate, in basepairs
        sequences <- msprime$simulate(length=as.integer(num_basepairs),
                                      migration_matrix=mig_matrix,
                                      mutation_rate=params@values[1][[i]]$mutation_rate,
                                      population_configurations=c(pop_meta, pop_local),
                                      demographic_events=c(split_event))
        
        # divide ts.diversity() by sequence length to get average number of pairwise difference __per basepair__
        bp_diversity <- sequences$diversity()/num_basepairs
      }
      g_diversities <- c(g_diversities,bp_diversity)
    }
    role@modelRuns[[i]]@timeseries[[niter_timestep]]@localComm@gdiversitiesSp <- g_diversities
  }
}

# population size -> ancestry -> ancestry + mutation rate -> gene sim regions 
# demography class - simple bc single population
# add population 
# record abundance every so often - take harmonic mean - plug in as initial size
# sequence length - 500 bp 
# samples - num individuals per population 
# for mess its like 10 

# same thing for 
# in mess, keep track of number of migration eevents, scale to sim, average migration rate
# in demography if multiple populations

# for every timesetep backward
# for eveyr species

# TODO
# per species simulate ancestry and get tree
# simulate mutations on tree and get tree sequence
# call diversity on tree sequence 
# outcome is diversity of each species in local community 
# create a demography of 1 population of J individuals 

# end up with vector of 0s and 1s, saying if same as reference or different from reference
# save