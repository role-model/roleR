library(roleR)
library(ggplot2)
library(dplyr)
set.seed(1988) #jbp
### Setting up a roleExperiment

#### Replicates of the same parameter settings

for(i in 1:100){
    print(i)
    params1 <- untbParams(individuals_local = 100, individuals_meta = 1000, 
                          species_meta = 50, 
                          speciation = 0.2, 
                          dispersal_prob = 0.1, init_type = 'oceanic_island',
                          niter = 10000, niterTimestep = 100) 
    
    paramsList1 <- list(a= params1, b= params1, c=params1, d= params1, e=params1)
    
    trial1 <- roleExperiment(paramsList1)
    trial1 <- runRole(trial1)
}