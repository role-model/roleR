
#' @title An S4 class to specify one RoLE model - or 1 run within roleSim
#'
#' @slot timeseries a list of roleData objects
#' @slot paramValues a named list of numeric vectors containing params
#'
#' @export

setClass('roleModel',
         slots = c(timeseries = 'list', 
                   paramValues = 'list'))

# constructor for roleModel
roleModel <- function(timeseries,paramValues,niter)
{
  new("roleParams", timeseries=timeseries, paramValues=paramValues)
}

# "as" functions don't work with C++ pointers, so created From and To functions
roleModelFromCpp <- function(modelCpp) {
  
  ts <- list() 
  
  for(i in 1:length(modelCpp$timeseries))
  {
    ts <- append(ts,roleDataFromCpp(modelCpp$timeseries[[i]]))
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

# create a model from scratch with default params, exclusively for testing purposes
dummyModel <- function(R=TRUE, run=FALSE){
  
  params <- roleParams(nrun=1,niter=1000,niterTimestep=100,defaults=TRUE)
  cparams <- stretchAndSampleParams(params)
  parlist <- cparams@values[[1]]
  model <- initSim(parlist,type="bridge_island",niter=1000)

  model$print <- FALSE
  model$local$print <- FALSE

  if(run){
    iterSim(model,params@niter,params@niterTimestep,FALSE)
    model
    model$
    model$timeseries[[1]]
    model$local
    model$timeseries
    test <- model$timeseries
    
  }
  if(R){
    model <- roleModelFromCpp(model)
  }
  return(model)
}
