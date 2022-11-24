model@modelSteps[[1]]@metaComm@spAbund
model@modelSteps[[1]]@localComm@spAbund

model@params

tre <- as(model@modelSteps[[2]]@phylo, 'phylo')
plot(tre)
model@modelSteps[[2]]@phylo@tipNames

roledata <- model@modelSteps[[2]]


# should the trimming happen each time step, each timestep 
# should have a phylo with all the alive local species?
# exclude species in meta that are not alive in lcoal

getLocalPhylo <- function(roledata){
    species_meta <- length(roledata@metaComm@spAbund)
    sp_meta_names <- paste0("t",1:species_meta)
    sp_meta_names <- sp_meta_names[(!roledata@localComm@spAbund > 0)[1:50]]
    phy <- as(roledata@phylo, 'phylo')
    phy <- drop.tip(phy, sp_meta_names)
    return(phy)
}

tre$tip.label
species_meta <- 50
sp_meta_names <- paste0("t",1:species_meta)
sp_meta_names <- sp_meta_names[(!model@modelSteps[[2]]@localComm@spAbund > 0)[1:50]]

test <- drop.tip(tre, sp_meta_names)
plot(test)
(!model@modelSteps[[2]]@localComm@spAbund > 0)[1:50]
plot(test)
drop.tip(tre, paste0("t",1:50), trim.internal = TRUE, subtree = FALSE,
         root.edge = 0, rooted = is.rooted(phy), collapse.singles = TRUE,
         interactive = FALSE)