
roleModelFromCpp <- function(modelCpp) {
  localComm <- localComm(modelCpp$local$abundance_indv,modelCpp$local$species_ids,
                         modelCpp$local$traits,modelCpp$local$abundance_sp,
                         modelCpp$local$traits_sp,modelCpp$local$pi)
  return(new('roleModel',
             localComm = localComm))
}

#' @title An S4 class to specify the entire RoLE model
#'
#' @slot timeseries an object of class \code{list} of roleData objects
#' @slot params a named list of numeric vectors containing params
#' 
#' @export

setClass('roleModel',
         slots = c(timeseries = 'list',params = 'list'))

#' @title Create a RoLE model object
#'
#' @param modelCpp a roleModelCpp object 
#'
#' @return an object of class \code{roleModel}
#'
#' @export

roleModel <- function(modelCpp) {
  localComm <- localComm(modelCpp$local$abundance_indv,modelCpp$local$species_ids,
                         modelCpp$local$traits,modelCpp$local$abundance_sp,
                         modelCpp$local$traits_sp,modelCpp$local$pi)
  return(new('roleModel',
             localComm = localComm))
}