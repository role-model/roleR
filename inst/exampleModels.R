# exampleModels contains examples of created sim models
library(roleR)
sim <- initSim()
sim2 <- sim$copy()
sim2$death()


sim$local$abundance_indv
sim2$local$abundance_indv

sim$print = FALSE
sim$local$print = FALSE
iterSim(sim,nsteps=100,50,print=FALSE)


sim$local$species_ids
sim$local$abundance_sp
View(sim$local$traitdiffs)

# traitdiffs and abundance_sp not updating properly 

sim <- iterSim(sim,1,true)
sim$local$species_ids

20%%20
