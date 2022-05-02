#' @title An S4 class to role model data for timeseries
#' @slot localComm an object of class \code{localComm}
#' @slot metaComm an object of class \code{metaComm}
#' @slot phylo an object of class \code{rolePhylo}
#' @slot stats an object of class \code{data-frame} - 1st col contains param names,
#' 2nd col contains numeric entropies used to compute hill statistics, 
#' and 3rd col contains hill statistic values
#' 
#' @export

setClass('roleData',
         slots = c(localComm = 'localComm',
                   metaComm = 'metaComm',
                   phylo = 'rolePhylo',
                   stats = 'data.frame',
                   iterNum = 'integer'))

# constructor
roleData <- function(localComm,metaComm,phylo,iterNum) {
  return(new('roleData',
             localComm = localComm,
             metaComm = metaComm,
             phylo = phylo,
             iterNum = iterNum))
}

roleDataFromCpp <- function(data) {
  local <- localComm(data$local$abundance_indv,data$local$species_ids,
                     data$local$traits,data$local$abundance_sp,
                     data$local$traits_sp,data$local$pi_sp)
  
  meta <- metaComm(data$meta$abundance,data$meta$traits)
  phylo <- rolePhyloFromCpp(data$phylo)
  stats <- data$stats 
  iterNum <- data$iter_num
  
  return(roleData(local,meta,phylo,iterNum))
}

roleDataToCpp <- function(data) {
  local <- localComm(data$local$abundance_indv,data$local$species_ids,
                     data$local$traits,data$local$abundance_sp,
                     data$local$traits_sp,data$local$pi_sp)
  
  meta <- metaComm(data$meta$abundance,data$meta$traits)
  phylo <- rolePhyloFromCpp(data$phylo)
  stats <- data$stats 
  iterNum <- data$iter_num
  
  return(roleData(local,meta,phylo,iterNum))
}