#' @title A collection of data representing one state of a roleModel.
#' 
#' @description An S4 class to represent a roleModel state using a `localComm`,`metaComm`, & `rolephylo`.
#' 
#' @slot localComm A `localComm` object containing the local community.
#' @slot metaComm A `metaComm` object containing the meta community.
#' @slot phylo A `rolePhylo` object containing the model phylogeny. 
#' 
#' @details access the `localComm`,`metaComm`, and `phylo` slots using `@`, 
#' i.e. `data@localComm` for a `roleData` object called data 
#' @rdname roleData
#' @include rolePhylo.R
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