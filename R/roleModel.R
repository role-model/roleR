
#' @title An S4 class to specify one RoLE model - or 1 run within roleSim
#'
#' @slot timeseries a list of roleData objects
#' @slot paramValues a named list of numeric vectors containing params
#'
#' @export

setClass('roleModel',
         slots = c(timeseries = 'list', 
                   params = 'roleParams',
                   niter = 'numeric'))

# constructor for roleModel
roleModel <- function(timeseries,params,niter)
{
  new("roleModel", timeseries=timeseries, paramValues=paramValues,niter=niter)
}

# "as" functions don't work with C++ pointers, so created From and To functions
roleModelFromCpp <- function(modelCpp) {
  
  ts <- list() 

  if(length(modelCpp$timeseries) > 0){
  for(i in 1:length(modelCpp$timeseries))
  {
    ts <- append(ts,roleDataFromCpp(modelCpp$timeseries[[i]]))
  }
  }
  
  return(new('roleModel',
             timeseries=ts,
             paramValues = modelCpp$params))
}

roleModelToCpp <- function(model) {
  localComm <- localComm(modelCpp$local$abundance_indv,modelCpp$local$species_ids,
                         modelCpp$local$traits,modelCpp$local$abundance_sp,
                         modelCpp$local$traits_sp,modelCpp$local$pi)
  return(new('roleModel',
             localComm = localComm))
}

#niters = 100
# create a model from scratch with default params, exclusively for testing purposes
# dummyModel example runs:
# role <- dummyModel(run=TRUE,niters=1000,return_experiment=TRUE)
dummyModel <- function(R=TRUE, run=FALSE,fill_ts=FALSE,niters=100,return_experiment=FALSE,print=FALSE,no_speciation=F,no_dispersal=F){
  
  args <- roleArguments(nruns=1,niter=niters,niterTimestep=niters/10,defaults=TRUE)
  #args <- roleArguments(nruns=1,niter=100,niterTimestep=10,defaults=TRUE)
  
  if(no_dispersal){
    args@params[[1]]@dispersal_prob <- 0
  }
  cargs <- stretchAndSampleParams(args)
  parlist <- cargs@params[[1]]
  model <- initModel(parlist,type="bridge_island",niter=niters)
  
  if(no_speciation){
    model$params$speciation_local <- 0
    model$params$speciation_meta <- 0
  }

  if(fill_ts){
    data_copy <- model$copyData() 
  }
  
  if(print == FALSE){
    model$print <- FALSE
    model$local$print <- FALSE
  }

  if(run){
    model <- iterSim(model,args@niter,args@niterTimestep,print)
  }
  
  if(fill_ts){
    ts <- vector("list", niters/(niters/10))
    for(i in 1:length(ts)){
      ts[[i]] <- data_copy
    }
    model$timeseries <- ts
  }
  
  if(R){
    model <- roleModelFromCpp(model)
  }
  
  if(return_experiment){
    runs = list(model) 
    return(new("roleExperiment", modelRuns=runs, arguments=args))
  }
  else{
    return(model)
  }
}
