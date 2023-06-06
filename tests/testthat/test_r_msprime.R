# reticulate::conda_list()
condaenv_name <- "r-reticulate"
reticulate::use_condaenv(condaenv_name)
#py_list_packages()
msprime <- reticulate::import('msprime')
newick <- reticulate::import('newick')
collections <- reticulate::import('collections')



library(roleR)

p <- roleParams(dispersal_prob = 0.5, alpha = 5000, individuals_local = 1000, 
                niter = 15000, niterTimestep = 5000)

m <- runRole(roleModel(p))


