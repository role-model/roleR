# exampleModels contains examples of created sim models
library(roleR)
library(Rcpp)

sim <- initSim()
sim
sim$print = FALSE
sim$local$print = FALSE
iterSim(sim,10000,print=FALSE)


sim$local$species_ids
sim$local$abundance_sp
View(sim$local$traitdiffs)

# traitdiffs and abundance_sp not updating properly 

sim <- iterSim(sim,1,true)
sim$local$species_ids

# install.packages("devtools")
devtools::install_github("rdinnager/slimr")