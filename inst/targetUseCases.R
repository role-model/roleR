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

p <- roleParams(nsim=100,niter=1000)


setParam(p,param_name="speciation_meta",value=0.15,sim_num=1) # sets the speciation_meta param to 0.15 for the sim run #1

v <- runif(p@niter,min=0,max=0.3) # sample values from a unif distribution
setParam(p,param_name="speciation_meta",value=v,sim_num=1) # set values for sim run #1 
setParam(p,param_name="speciation_meta",value=v) # set values for all sim runs 

r@params[[1]]$speciation_meta
r <- addDefaultParams(r) 

saveRDS(r, file="x_out.Rda")
x_out2 <- readRDS(file="x_out.Rda")

setParam(r,"param")
setParam(r,"param")
rep(v,100)
v <- 0

setParam(r,"sigma_bm")
param_name <- "speciation_local"
r@params[[1]][,param_name] <- rep(1:1000)
r@params[[1]][,param_name]

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



