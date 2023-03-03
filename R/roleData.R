#' @title A collection of data representing one state of a roleModel
#' 
#' @description An S4 class to represent a roleModel state using a `localComm`,`metaComm`, & `phylo`
#' 
#' @slot localComm a `localComm` object containing the local community
#' @slot metaComm a `metaComm` object containing the meta community
#' @slot phylo a `rolePhylo` object containing the model phylogeny 
#' 
#' @rdname roleData
#' @export

setClass('roleData',
         slots = c(localComm = 'localComm',
                   metaComm = 'metaComm',
                   phylo = 'rolePhylo'))

# constructor
#' @rdname roleData
#' @export

roleData <- function(localComm,metaComm,phylo) {
  return(new('roleData',
             localComm = localComm,
             metaComm = metaComm,
             phylo = phylo))
}

roleDataFromCpp <- function(data) {
  local <- localComm(data$local$abundance_indv,data$local$species_ids,
                     data$local$traits,data$local$abundance_sp,
                     data$local$traits_sp,data$local$pi_sp)
  
  meta <- metaComm(data$meta$abundance,data$meta$traits)
  phylo <- rolePhyloFromCpp(data$phylo)

  return(roleData(local,meta,phylo))
}

roleDataToCpp <- function(data) {
  local <- localComm(data$local$abundance_indv,data$local$species_ids,
                     data$local$traits,data$local$abundance_sp,
                     data$local$traits_sp,data$local$pi_sp)
  
  meta <- metaComm(data$meta$abundance,data$meta$traits)
  phylo <- rolePhyloFromCpp(data$phylo)
  
  return(roleData(local,meta,phylo))
}