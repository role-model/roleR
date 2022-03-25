library(roleR)

# --------------
# create a simulation with default parameters and metadata, get and plot various values
# --------------

params <- roleParams(nruns=10,niter=1000)
sim <- roleSim(params)

# look at the total number of species that have existed in the model
state <- getStateData(sim,iter=1,iterAsPercentile=TRUE) # get the end state (100th percentile of iters) of the first run 
state <- getStateData(sim,runNum=2,iter=0.5,iterAsPercentile=TRUE) # get the mid state (50th percentile of iters) of the second run 

# add Hill numbers to all runs in the sim 
addHillStats(sim)

# get the trait diversity value
state@stats["trait_diversity"]

# get a list of species abundances over time
popseries <- getTimeseries(sim,type="model_value",runNum=1,valueName="abundanceSp")

# get the model phylogeny
library(ape)
endPhylo <- getEndState(sim@modelRuns[1])@phylo
plot(as(endPhylo, ape::rphylo)) # plot the ape phylogeny

# --------------
# create and modify parameters, priors, and hyperpriors
# --------------

# create a default set of params
p <- roleParams(nsim=3,niter=1000)

# create a default set of params with some individual param modifications
# allow params here to be a list 
p <- roleParams(nsim=3,niter=1000,speciation_meta=0.15)

# set a param 
p <- setParam(p,paramName="speciation_meta",value=0.15,runNum=1) # sets the speciation_meta param to 0.15 for the sim run #1
p <- setParam(p,param_name="speciation_meta",value=0.15,runNum=c(1,2,3)) # sets the speciation_meta param to 0.15 for the sim runs #1, #2, and #3
p <- setParam(p,param_name="speciation_meta",value=runif(1000,0,0.3),runNum=1) # sets the speciation_meta param to be drawn from a unif for the sim run #1
p <- setParam(p,param_name="speciation_meta",value=list(0.15,0.2,0.3),runNum=c(1,2,3)) # sets the speciation_meta param to 0.15,0.2,0.3 for runs 1 through 3 

# set a prior
unif_prior <- function(n,min=0,max=0.5){ # create a uniform prior sampling between 0 and 0.5 
  return(runif(n=n,min=min,max=max))}
p <- setPrior(p, param_name="speciation_local",priorFunc=unif_prior) # set the prior for the first two runs
p <- setPrior(p, param_name="speciation_local",priorFunc=list(unif_prior,unif_prior,unif_prior), runNum = c(1,2)) # set the prior for all runs

# set an iterfunc
# may want to call this setHyperPrior if that's clearer
p <- setIterFunc()

# --------------
# OLD BELOW

# set different time varying parameters for each sim_num 
# set different static params for each sim num 
# function that takes an iteration and returns a value
# a single value
# a single function
# or a list of values 
# or a list of functions
p <- 

v <- runif(p@niter,min=0,max=0.3) # sample values from a unif distribution
setParam(p,param_name="speciation_meta",value=v,sim_num=1) # set values for sim run #1 
v <- runif(p@niter,min=0,max=0.3) # sample values from a unif distribution
setParam(p,param_name="speciation_meta",value=v) # OR set values for all sim runs 
setPrior(p, param_name="speciation_local",priorFunc = runif, value=v) # OR just set the prior

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


# --------------
# create a simulation specifying parameter priors 
# --------------

