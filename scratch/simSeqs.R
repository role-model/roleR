
simSeqs <- function(model){
    
    # convert to ape phy
    phy <- as(model@modelSteps[[1]]@phylo, 'phylo')
    #phy <- ape::rphylo(10, p@speciation_meta, p@extinction_meta,T0=0)
    #plot(phy)
    #plot(phy,use.edge.length = TRUE)
    #add.scale.bar()
    #axisPhylo()
    #edgelabels(phy$edge.length,cex=0.5)
    #phy$edge.length <- rep(0,length(phy$edge.length))
    
    # create demography object
    # initial size should be the initial size of each pop
    d <- msprime$Demography$from_species_tree(ape::write.tree(phy), 
                                              time_units = 'gen', 
                                              initial_size = model@params@individuals_meta / model@params@species_meta, 
                                              growth_rate = 0)
    # create sp abundance matrix
    abundance_ts_mat <- matrix(data=0,nrow=(model@params@niter/model@params@niterTimestep) + 1,ncol=10000)
    for(m in 1:length(model@modelSteps)){
        abundance_ts_mat[m,] <- model@modelSteps[[m]]@localComm@spAbund
    }
    
    sp_indices <- which(colSums(abundance_ts_mat) > 0)
    
    #sp_index = 1
    for(sp_index in sp_indices){
        
        # get an abundance vector 
        sp_abundance_vect <- abundance_ts_mat[,sp_index]
        
        # SCENARIOS
        # 1. speciation from meta - a split of the meta pop into a new local and meta pop
        # Check if species came from a meta ancestor (ancestor has a tipname < n meta species)
        
        # 2. speciation from local - a split of the local pop into two new local pops
        # Check if species came from a local ancestor (ancestor has a tipname > n meta species)
        
        # 3. immigration from meta - a split of the meta pop into a new local and meta pop
        # Check if species came from meta (sp index < n meta species)
        
        # 4. birth from local - don't care about this because it's not a new pop
        
        # things can immigrate, go extinct, then immigrate again
        
        # emerges holds all new emergences of the species in the local
        # or times when the abund went from 0 to >0 
        emerges <- c()
        
        # get emergence at timestep 1
        if(sp_abundance_vect[1] > 0){
            emerges <- c(emerges,1)
        }
        # get all other emergences
        for(ts_i in 2:length(sp_abundance_vect)){
            if(sp_abundance_vect[ts_i-1] == 0 & sp_abundance_vect[ts_i] > 0){
                emerges <- c(emerges,ts_i)
            }
        }
        
        # for each emergence
        for(e in emerges){
            
            # immigration from meta
            # if sp index < n species meta
            if(sp_index < model@params@species_meta){
                
                pop_name <- paste0('t', as.character(sp_index))
                pop_name_l <- paste0(pop_name, '_l')
                pop_name_m <- paste0(pop_name, '_m')
                d$add_population(name = pop_name_l,initial_size=100)
                d$add_population(name = pop_name_m, initial_size = 100)
                d$add_population_split(time = e, derived = c(pop_name_l, pop_name_m), 
                                       ancestral = pop_name)
            }
            
            # speciation from meta 
            # if parent node index < n species meta
            else if(which(model@modelSteps[[11]]@phylo@e[,2] == sp_index) < model@params@species_meta)
            {
                pop_name_l <- paste0('t', as.character(sp_index))
                pop_name_m <- paste0('t',which(model@modelSteps[[11]]@phylo@e[,2]))
                
                d$add_population(name = pop_name_l,initial_size=1)
                d$add_population(name = pop_name_m, initial_size = 100)
                d$add_population_split(time = e, derived = c(pop_name_l, pop_name_m), 
                                       ancestral = pop_name_m)
            }
            
            # speciation from local
            else{
                pop_name_child <- paste0('t', as.character(sp_index))
                pop_name_anc <- paste0('t',which(model@modelSteps[[e]]@phylo@e[,2]))
                
                d$add_population(name = pop_name_l,initial_size=1)
                d$add_population(name = pop_name_m, initial_size = 100)
                d$add_population_split(time = e, derived = c(pop_name_l, pop_name_m), 
                                       ancestral = pop_name_m)
            }
        }
        
        ts <- msprime$sim_ancestry(reticulate::dict(list(A_meta = 3, A_local = 4, B = 2, C = 1)), demography = d, ploidy = 1)
        #first_nonzero <- min(which(sp_abundance_vect > 0))
    }
        #model@modelSteps[[i]]@timeseries[[niter_timestep]]@localComm@gdiversitiesSp <- g_diversities
}


#msprime <- import('msprime')
simSeqs2 <- function(model){
    library(plyr)
    
    # diem_scalar <- 1 
    
    # convert to ape phy
    phy <- as(model@modelSteps[[length(model@modelSteps)]]@phylo, 'phylo')
    
    
    # create demography object
    # error - all leaf pops must be at time 0 
    # initial size should be the initial size of each pop
    d <- msprime$Demography$from_species_tree(ape::write.tree(phy), 
                                              time_units = 'gen', 
                                              initial_size = 1000, growth_rate = 0)
    
    # create sp abundance matrix
    abundance_ts_mat <- matrix(data=0,nrow=(model@params@niter/model@params@niterTimestep) + 1,ncol=10000)
    for(m in 1:length(model@modelSteps)){
        abundance_ts_mat[m,] <- model@modelSteps[[m]]@localComm@spAbund
    }
    
    # create sp last origin step matrix
    origin_step_ts_mat <- matrix(data=0,nrow=(model@params@niter/model@params@niterTimestep) + 1,ncol=10000)
    for(m in 2:length(model@modelSteps)){
        origin_step_ts_mat[m,] <- model@modelSteps[[m]]@localComm@spLastOriginStep
    }
    
    # create sp extinction step matrix
    ext_step_ts_mat <- matrix(data=0,nrow=(model@params@niter/model@params@niterTimestep) + 1,ncol=10000)
    for(m in 2:length(model@modelSteps)){
        ext_step_ts_mat[m,] <- model@modelSteps[[m]]@localComm@spExtinctionStep
    }
    
    # create sp harm mean matrix
    harm_ts_mat <- matrix(data=0,nrow=(model@params@niter/model@params@niterTimestep) + 1,ncol=10000)
    for(m in 2:length(model@modelSteps)){
        harm_ts_mat[m,] <- model@modelSteps[[m]]@localComm@spAbundHarmMean
    }
    
    # get all species that existed in the local at some point 
    sp_indices <- which(colSums(abundance_ts_mat) > 0)
    
    # create function to be used in finding next extinction after an origin
    getNextPositiveIndexAfterIndex <- function(v,index){
        for(e in v[index:length(v)]){
            if(e != 0){
                return(e)
            }
        }
        return(-1)
    }
    
    # for each species
    for(sp_index in sp_indices)
    {
        # get values for this species 
        sp_abundance_vect <- abundance_ts_mat[,sp_index]
        origin_vect <- origin_step_ts_mat[,sp_index]
        ext_vect <- ext_step_ts_mat[,sp_index]
        harm_vect <- harm_ts_mat[,sp_index]
        
        # an origin is a timestep when the abundance changed from 0 to 1 
        # for each unique origin
        uniq_origins <- unique(origin_vect)
        for(o in uniq_origins){
            
            # get the following extinction
            ext_index <- getNextPositiveIndexAfterIndex(ext_vect,o)
            
            # get the name of the pop tip 
            # I believe species indices follow phylo tip names? or should do so? 
            pop_name <- paste0('t', as.character(sp_index))
            # create a new name for the local and meta
            pop_name_l <- paste0(pop_name, '_l')
            pop_name_m <- paste0(pop_name, '_m')
            d$add_population(name = pop_name_l,initial_size=1)
            approx_meta_abund <- model@modelSteps[round_any(o, model@params@niterTimestep)]
            d$add_population(name = pop_name_m,initial_size=approx_meta_abund)
            d$add_population_split(time = o, derived = c(pop_name_l, pop_name_m), 
                                   ancestral = pop_name)
            # add the extinction time
        }
    }
}