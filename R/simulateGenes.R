# install miniconda, r-reticulate, and msprime
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

roleulateSpeciesGeneticDiv <- function(role){
  
  library(reticulate)
  msprime <- import("msprime")
  
  # get a model from role experiment
  params <- role@params
  model <- role@modelRuns[[1]]
  niter <- role@niter
  niter_timestep <- role@niter_timestep
  num_basepairs <- params@values[1][[1]]$num_basepairs

  # get the model timeseries for species abundances
  #abundance_ts <- getTimeseries(run,type="model_value",valueName="abundanceSp")
  abundance_ts <- getTimeseries(run,type="model_value",valueName="abundanceSp")
  
  # convert to a matrix 
  # delete cols where every value is 0 
  
  genetic_diversities <- numeric()
  
  # for each species
  for(sp_index in 1:ncol(abundance_ts))
  {
    # get vector of abundances over time for that species
    abundance_vect <- abundance_ts[,sp_index]
      
    # get index of first nonzero abundance - this is the timestep at which the species started
    start_timestep <- min(which(abundance_vect != 0))
    end_timestep <- max(which(abundance_vect != 0))
    
    n_steps_lived <- end_timestep - start_timestep
      
    local_pop_start <- abundance_vect[start_timestep]
    local_pop_end <- abundance_vect[end_timestep]
    
    # if the species originated in local and thus did not come from meta
    if(sp_index > model@timeseries[[1]]@metaComm@abundance)
    {
      # sample size is species abundance at present
      # Ne is effective populatin size - harmonic mean of all abundances of all timesteps
      # Ne is alpha times harmonic mean 
      # alpha gets estimated pretty accurately
      # below no meta
      ts <- msprime.simulate(sample_size=10,
                             Ne=1e5,
                             length=num_basepairs,
                             mutation_rate=1e-7)
      bp_diversity <- ts.diversity()/num_basepairs
    }
    # otherwise the species immigrated from meta
    else
    {
      meta_pop_startend <- model@timeseries[[1]]@metaComm@abundance[sp_index] * 100 # * individuals meta
      
      # create the meta community
      pop_meta <- msprime$PopulationConfiguration(sample_size=meta_pop_startend,
                                                  initial_size=meta_pop_startend)
      pop_local <- msprime$PopulationConfiguration(sample_size=local_pop_end,
                                                   initial_size=local_pop_start)
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
                                           proportion=1)
      
    
      # Length of sequence to simulate, in basepairs
      sequences <- msprime$simulate(length=as.integer(num_basepairs),
                             migration_matrix=mig_matrix,
                             mutation_rate=params@values[1][[1]]$mutation_rate,
                             population_configurations=c(pop_meta, pop_local),
                             demographic_events=c(split_event))
      
      # divide ts.diversity() by sequence length to get average number of pairwise difference __per basepair__
      bp_diversity <- ts.diversity()/length
    }
  }
  genetic_diversities <- c(genetic_diversities,bp_diversity)
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