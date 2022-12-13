#' @title using msprime simulation, add genetic diversities to a roleExperiment or roleModel for each species 
#'
#' @param x a roleModel or roleExperiment object to get a stat from 

#' @return a dataframe containing the timeseries 
#'
#' @export

#alpha = 0.5
simulateSpeciesGenDiv <- function(model){
    
    # import msprime into R 
    library(reticulate)
    msprime <- import("msprime")
    
    # init spGenDiv with zeroes
    for(s in 1:length(model@modelSteps)){
        model@modelSteps[[s]]@localComm@spGenDiv <- rep(0,10000)
    }
    
    # get params from model
    diem_scalar <- 100
    niter <- model@params@niter
    niter_timestep <- model@params@niterTimestep
    indv_local <- model@params@individuals_local[1]
    alpha <- model@params@alpha # alpha is not a slot yet 
    num_basepairs <- model@params@num_basepairs
    J <- model@params@individuals_local[1]
    
    # create sp abundance matrix
    abundance_ts_mat <- matrix(data=0,nrow=(niter/niter_timestep) + 1,ncol=10000)
    for(m in 1:length(model@modelSteps)){
        abundance_ts_mat[m,] <- model@modelSteps[[m]]@localComm@spAbund
    }
    
    # delete cols where every value is 0 - don't understand why this is again, ask 
    # resolves problem of resolution where a species could exist then go extinct within one interval
    
    # get species indices to loop through
    sp_indices <- which(colSums(abundance_ts_mat) > 0)
    
    # initialize g_diversities
    #g_diversities <- numeric(ncol(abundance_ts))
    # place NAs for species where no individuals are recorded to have existed 
    #g_diversities[-sp_indices] <- NA
    
    # sp_index = 1
    # for each species (each column in the abundance timeseries matrix) 
    for(sp_index in sp_indices)
    {
        sp_abundance_vect <- abundance_ts_mat[,sp_index]
        
        # diem size - scalar for difference between observed and actual counts 
        sp_abundance_vect <- sp_abundance_vect * diem_scalar
        
        # for every timestep pair i.e. a step and that step + 1
        # t = 1
        for(t in 1:(length(sp_abundance_vect)-1)){
            
            if(sp_abundance_vect[t] > 0 & sp_abundance_vect[t+1] > 0){ #if(sp_abundance_vect[t] > 0 & sp_abundance_vect[t+1] > 0){
                
                # get the start and end timesteps to simulate
                start_timestep <- t
                end_timestep <- t + 1
                
                # start - end * iters per timestep
                # 1 generation = J/2 timesteps 
                n_gens_lived <- ((end_timestep - start_timestep) * niter_timestep) / (J/2)
                local_pop_start <- sp_abundance_vect[start_timestep]
                local_pop_end <- sp_abundance_vect[end_timestep]
                
                start_end_vect <- sp_abundance_vect[start_timestep:end_timestep]
    
                # if the species originated in local and thus did not come from meta
                if(sp_index > length(model@modelSteps[[1]]@metaComm@spAbund) | sp_index == 0)
                {
                    # sample size is species abundance at present
                    # Ne is effective populatin size - harmonic mean of all abundances of all timesteps
                    # Ne is alpha times harmonic mean 
                    # alpha gets estimated pretty accurately
                    # below no meta
                    
                    # num_basepairs = 10
                    ts <- msprime$simulate(sample_size=as.integer(10),
                                           Ne=length(start_end_vect) / sum(1/start_end_vect),
                                           length=as.integer(num_basepairs),
                                           mutation_rate=as.integer(1e-7))
                    bp_diversity <- ts$diversity()/num_basepairs

                    # at the timeseries, for the given species, set bp_diversity
                    model@modelSteps[[t+1]]@localComm@spGenDiv[sp_index] <- bp_diversity
                }
                
                # otherwise the species immigrated from meta
                else
                {
                    meta_pop_startend <- model@modelSteps[[1]]@metaComm@spAbund[sp_index] * 100 # * individuals meta
      
                    # create the meta community
                    pop_meta <- msprime$PopulationConfiguration(sample_size=as.integer(meta_pop_startend),
                                                                initial_size=as.integer(meta_pop_startend))
                    pop_local <- msprime$PopulationConfiguration(sample_size=as.integer(local_pop_end),
                                                                 initial_size=as.integer(local_pop_start))
                    
    
                    # specify migration from local to metacommunity (backwards in time)
                    #mig_matrix <- matrix(c(0,arguments@values[1][[1]]$dispersal_prob,0,0),nrow=2,ncol=2)
                    mig_matrix <- matrix(c(0,model@params@dispersal_prob[1],0,0),nrow=2,ncol=2) #todo fix
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
                                                  mutation_rate=model@params@mutation_rate[1], #todo fix
                                                  population_configurations=c(pop_meta, pop_local),
                                                  demographic_events=c(split_event))
                    
                    # divide ts.diversity() by sequence length to get average number of pairwise difference __per basepair__
                    bp_diversity <- sequences$diversity()/num_basepairs
                    # at the timeseries, for the given species, set bp_diversity
                    model@modelSteps[[t+1]]@localComm@spGenDiv[sp_index] <- bp_diversity
                }
                #g_diversities <- c(g_diversities,bp_diversity)
            }
        }
        #model@modelSteps[[i]]@timeseries[[niter_timestep]]@localComm@gdiversitiesSp <- g_diversities
    }
    return(model)
}