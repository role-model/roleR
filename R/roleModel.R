
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
    append(ts,roleDataFromCpp(modelCpp$timeseries(i)))
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
