library(roleR)

# --------------
# create a simulation with default parameters and metadata 
# --------------
sim <- roleSim() # same as sim <-roleSim(params = NULL, init = NULL, nstep = 100, series_timestep = NULL, nsim = 1, print = FALSE) 

# look at some values 
sim$local$Imax #the max number of individuals in the local community
plot(getApePhylo(sim)) # plot the ape phylogeny
plot(sim$local$abundance_sp) # plot species abundances 


# --------------
# create a simulation specifying distributions of parameters 
# --------------

params <- roleParameters()
params@timeseries$dispersal_prob <- runif(10000,min=0,max=1) # dispersal_prob is now randomly sampled from a uniform distribution

# timeseries is a sequence
# distribution is in parallel (10000 randomly initialized communitii) could be for one param 
# combinations replicate many or for different param combinations replicate once
# create a function for exponential increase of a parameter from 0 to 0.6
# parallel timeseries 

eq <- function(x){x*x}
curve(eq, from=0, to=1, , xlab="x", ylab="y") #plot it 
series <- eq(1:1000) / 1000000 * 0.6 
params@timeseries$speciation_local <- series
