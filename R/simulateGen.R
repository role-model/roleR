simulateGen <- function(model){
    
    # grab static params
    local_pop <- model@params@individuals_local[1]
    meta_pop <- model@params@individuals_meta[1]
    
    niter <- model@params@niter
    niter_timestep <- model@params@niterTimestep
    
    alpha <- model@params@alpha # alpha is not a slot yet 
    bp <- model@params@num_basepairs
    ploidy <- 2
    diem_scalar <- 100
    mutation_rate <- 0.01
    
    # create sp abundance matrix
    abundance_ts_mat <- matrix(data=0,nrow=(niter/niter_timestep) + 1,ncol=10000)
    for(m in 1:length(model@modelSteps)){
        abundance_ts_mat[m,] <- model@modelSteps[[m]]@localComm@spAbund
    }
    #abundance_ts_mat <- abundance_ts_mat[-1,]
    
    # delete cols where every value is 0 - don't understand why this is again, ask 
    # resolves problem of resolution where a species could exist then go extinct within one interval
    
    # get species indices to loop through
    sp_indices <- which(colSums(abundance_ts_mat) > 0)
    sp_index <- 1
    for(sp_index in sp_indices)
    {
        # diem size - scalar for difference between observed and actual counts 
        sp_abundance_vect <- abundance_ts_mat[,sp_index] * diem_scalar
        
        for(i in 2:length(model@modelSteps)){
            local_samp <- sp_abundance_vect[i]
            if(local_samp != 0 & local_samp > 1){
                imm_rate <- model@params@dispersal_prob[(i-1) * niter_timestep]
                gen <- simGen(local_pop,meta_pop,local_samp, mutation_rate,imm_rate,bp=bp,maxTime = 100, ploidy=ploidy)
            }
        }
    }
}