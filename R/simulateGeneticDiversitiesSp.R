
#' @title using msprime simulation, add genetic diversities to a roleExperiment or roleModel for each species 
#'
#' @param x a roleModel or roleExperiment object to get a stat from 
#' @param runNum an integer specifying the run number if a roleExperiment object is used
#' @param type a string specifying the type of values to get, either "summary_stats", or "model_values" 

#' @return a dataframe containing the timeseries 
#'
#' @export

simulateGeneticDiversitiesSp <- function(role){
  
  # NOTE - where does diemScalar get specified?
  
  # import msprime into R 
  # library(reticulate)
  msprime <- import("msprime")
  
  # get arguments from role experiment
  arguments <- role@arguments
  niter <- arguments@niter
  niter_timestep <- arguments@niterTimestep
  indv_local <- arguments@values[[1]][["individuals_local"]]
  alpha <- arguments@values[["alpha"]]
  
  # for every model run in the role experiment
  for(i in 1:length(role@modelRuns))
  {
    # i = 1
    # get model 
    model <- role@modelRuns[[i]]
    # get num basepairs
    num_basepairs <- arguments@values[[i]]$num_basepairs
    
    # get the model timeseries of species abundances
    ts <- getTimeseries(model,type="model_values")
    abundance_ts <- ts$abundanceSp
    
    # convert to a matrix 
    abundance_ts <- t(matrix(unlist(abundance_ts),ncol=niter_timestep))
    
    # delete cols where every value is 0 - don't understand why this is again, ask 
    # resolves problem of resolution where a species could exist then go extinct within one interval
    # abundance_ts <- abundance_ts[, colSums(abundance_ts) > 0]
    
    # get species indices to loop through
    sp_indices <- which(colSums(abundance_ts) > 0)
  
    # initialize g_diversities
    #g_diversities <- numeric(ncol(abundance_ts))
    # place NAs for species where no individuals are recorded to have existed 
    #g_diversities[-sp_indices] <- NA
    
    # for each species (each column in the abundance timeseries matrix) 
    for(sp_index in sp_indices)
    {
      # get index of first nonzero abundance - this is the timestep at which the species emerged
      speciation_timestep <- min(which(abundance_vect != 0))
      # likewise find the extinction
      extinction_timestep <- max(which(abundance_vect != 0))

      # diem size - scalar for different between observed and actual counts 
      abundance_vect <- abundance_ts[,sp_index] * diem_scalar
      
      # for every timestep pair i.e. a step and that step + 1
      for(t in 1:length(abundance_vect)){
        
        if(t >= speciation_timestep & t < extinction_timestep){
          # get the start and end timesteps to simulate
          start_timestep <- abundance_vect[t]
          end_timestep <- abundance_vect[t+1]
          
          # start - end * iters per timestep
          # 1 generation = J/2 timesteps 
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
            
            # at the timeseries, for the given species, set bp_diversity
            model@timeseries[[t+1]]@localComm@gDiversitiesSp[sp_index] <- bp_diversity
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
            mig_matrix <- matrix(c(0,arguments@values[1][[1]]$dispersal_prob,0,0),nrow=2,ncol=2)
            
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
            
            # arguments@values[1][[1]]$mutation_rate <- 0.1
            # Length of sequence to simulate, in basepairs
            sequences <- msprime$simulate(length=as.integer(num_basepairs),
                                          migration_matrix=mig_matrix,
                                          mutation_rate=arguments@values[1][[i]]$mutation_rate,
                                          population_configurations=c(pop_meta, pop_local),
                                          demographic_events=c(split_event))
            
            # divide ts.diversity() by sequence length to get average number of pairwise difference __per basepair__
            bp_diversity <- sequences$diversity()/num_basepairs
            # at the timeseries, for the given species, set bp_diversity
            model@timeseries[[t+1]]@localComm@gDiversitiesSp[sp_index] <- bp_diversity
          }
          #g_diversities <- c(g_diversities,bp_diversity)
        }
      }
      #role@modelRuns[[i]]@timeseries[[niter_timestep]]@localComm@gdiversitiesSp <- g_diversities
    }
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