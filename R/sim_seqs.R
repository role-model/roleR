#' Sim seqs
#' @description
#' Takes a roleModel and simulates seq data for it. 
#' 
#'
#' @param model a completed roleModel
#'
#' @return model with seqs added
#' @export
#'
#' @importFrom reticulate source_python
#' @importFrom ape write.tree
sim_seqs <- function(model) {
    sim_seqs_path <- system.file("python", "role_msprime.py", package = "roleR")
        
    reticulate::source_python(sim_seqs_path)
    
    #exp = as(model, "roleExperiment")
    
    ## the timesteps at which data were recorded
    ## TODO: Doing it this way because the as() above is broken atm
    iterations <- c(0, which(1:model@params@niter %% model@params@niterTimestep == 0))
    
    ## Global parameters
    J_m <- model@params@individuals_meta
    mu <- model@params@mutation_rate
    sequence_length <- model@params@num_basepairs
    
    ## For each snapshot
    for (idx in 1:length(model@modelSteps)){
        ## Fetch the data that we need to run the msprime simulation
        ## Current time in timesteps (necessary to get divergence counting from t_0)
        ## TODO: If as(model, "roleExperiment") works then you can do this instead.
        #curtime = exp@experimentMeta$iterations[[idx]]
        curtime = iterations[[idx]]
        
        ## Prior to iterFun it was done this way:
        ##J <- model@params@individuals_local[[idx]]
        J <- model@params@individuals_local(idx)
        alpha <- model@params@alpha(idx)
        
        # metaTree <- newick(model@modelSteps[[idx]]@phylo)
        metaTree <- as(model@modelSteps[[idx]]@phylo, 'phylo')
        metaTree <- ape::write.tree(metaTree, append = FALSE, digits = 30, 
                                    tree.names = FALSE)
        metaAbund <- model@modelSteps[[idx]]@metaComm@spAbund
        
        localAbund <- model@modelSteps[[idx]]@localComm@indSpecies
        spAbundHarmMean <- model@modelSteps[[idx]]@localComm@spAbundHarmMean
        localTDiv <- model@modelSteps[[idx]]@localComm@spLastOriginStep
        
        ## Returns a dataframe with rows for pi, TajD and genotypes
        res <- py_msprime_simulate(J_m, J, curtime, metaTree, metaAbund, localAbund,
                                   spAbundHarmMean, localTDiv, alpha, sequence_length, mu,
                                   verbose=FALSE)
        
        ## Update the model with the results
        model@modelSteps[[idx]]@localComm@spGenDiv = unlist(res["pi",])
        model@modelSteps[[idx]]@localComm@indSeqs = unlist(res["gtypes",])
        
    }
    return(model)
}
