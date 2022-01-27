# fit test
library(roleR)
library(randomForest)
library(hillR)

sim <- initSim()

fit <- function(input, simrun)
{
  simrun <- sim
  sim$params$values
  
  ts_values <- data.frame(matrix(ncol = 40, nrow = length(simrun$timeseries)))
  
  #each row is a timestep of a sim
  for(i in 1:length(simrun$timeseries))
  {
    s <- simrun$timeseries[i]
    # make row of params and summary stats
    row <- c(s$params$values[1],s$params$values[2])
    ts_values[i] <- row #add to data
  }
  
  rf <- randomForest(param ~ stat1 + stat2 + stat3 ... , data=ts_values, ..., subset, na.action=na.fail, test = input$param)
  rf <- out 
}

# get a timeseries of a parameter, summary statistic, population etc.
getParamTimeseries <- function(simrun, param)
{
  out <- numeric(1:length(simrun$timeseries))
  for(i in 1:length(simrun$timeseries))
  {
    out[i] <- simrun$timeseries[i]$params$values[param]
  }
  return(out)
}

# computes and adds Hill numbers to a roleModelCpp simulation object 
computeHillNumbers <- function(simrun, param)
{
  for(i in 1:length(simrun$timeseries)) # for each time step in the model
  {
    s <- simrun$timeseries[i]
    
    h <- hill_func(comm = t(data.frame(s$local$abundance_sp)), # calculate alpha_diversity_func
                      traits = t(data.frame(s$local$traits_sp)), q = 0)
    simrun$stats[alpha_diversity_func] <- h #assign the hill number
    
    h <- hill_func(comm = t(data.frame(s$local$abundance_sp)), # alpha_diversity_taxa
                   traits = t(data.frame(s$local$traits_sp)), q = 0)
    simrun$stats[alpha_diversity_taxa] <- h 
    
    h <- hill_func(comm = t(data.frame(s$local$abundance_sp)), # alpha_diversity_phy
                   traits = t(data.frame(s$local$traits_sp)), q = 0)
    simrun$stats[alpha_diversity_phylo] <- h 
  }
}

# equiation for abudnance, triats, phylogeny, treating genetic diversity as abundance 
# paramters for how many hill numbers log(2) scale between 0 and 5  

