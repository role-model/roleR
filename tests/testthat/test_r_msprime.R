cl <- reticulate::conda_list()
condaenv_name <- "r-reticulate"



if(!(condaenv_name %in% cl$name)) {
    reticulate::conda_create("r-reticulate")
    reticulate::conda_create(envname = 'foo', 
                             packages = c('msprime', 'newick', 'collections'))
}

reticulate::use_condaenv(condaenv_name)


library(roleR)


p <- roleParams(dispersal_prob = 0.5, alpha = 1000, individuals_local = 1000,
                speciation_local = 0, mutation_rate = 1e-06, num_basepairs = 500,
                species_meta = 100, individuals_meta = 200,
                niter = 3000, niterTimestep = 1000)



m <- roleModel(p)
plot(m@modelSteps[[1]]@metaComm@spAbund)

mrun <- runRole(m)




# 
# 
# lapply(mrun@modelSteps, function(x) x@localComm@spAbundHarmMean)
# 
# mrun@modelSteps[[4]]@localComm
# 
# 
# getSumStats(mrun)
