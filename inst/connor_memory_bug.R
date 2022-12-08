
# Connor's params threw same error on my computer
# testing with same params as Connor except with 10x less timesteps saved 
params <- untbParams(
    individuals_local = 900,
    individuals_meta = 9000,
    species_meta = 900,
    speciation = 0.3,
    dispersal_prob = 0.5,
    init_type = "oceanic_island",
    niter = 9000,
    niterTimestep = 100
)

exp <- roleModel(params)

# runs on my computer, suggesting the problem is in the saving/storage of timesteps
m <- iterModel(exp)

objSize <- function(obj){
    return(format(object.size(obj), units='auto'))
}

# all model sizes are very similar
objSize(m@modelSteps[[90]])
objSize(m@modelSteps[[45]])
objSize(m@modelSteps[[1]])

# local and meta are tiny 
objSize(m@modelSteps[[45]]@localComm)
objSize(m@modelSteps[[45]]@metaComm)

# phylo is MASSSSSSIVE 
objSize(m@modelSteps[[45]]@phylo)

big_phy <- m@modelSteps[[45]]@phylo


# e is big, l is big, alive big, tipnames big
objSize(big_phy@e)
objSize(big_phy@l)
objSize(big_phy@alive)
objSize(big_phy@tipNames)

# converting to an ape phy reduces size tremendously!!
# and plots just fine
ape_phy <- as(big_phy, 'phylo')
objSize(ape_phy)
plot(phy)

# raw e is 1,801,798 x 2 so unsurprising it is so big 
e <- m@modelSteps[[45]]@phylo@e
dim(e)

# converted ape phy is 4380 x 2 
dim(ape_phy$edge)

# number of species is 2191 
big_phy@n

# the very first phylo is also huge! meaning that speciation rate isn't causing this
objSize(m@modelSteps[[1]]@phylo)

# and the R phylo is as huge as the Cpp phylo!
objSize(exp@modelSteps[[1]]@phylo)
View(exp@modelSteps[[1]]@phylo@e)

# testing with small speciation rate but same timesteps
params <- untbParams(
    individuals_local = 900,
    individuals_meta = 9000,
    species_meta = 900,
    speciation = 0.001,
    dispersal_prob = 0.5,
    init_type = "oceanic_island",
    niter = 9000,
    niterTimestep = 10
)

exp <- roleModel(params)

# runs on my computer, suggesting the problem is in the saving/storage of timesteps
m <- iterModel(exp)

rm()
