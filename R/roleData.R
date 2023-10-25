#' @title A collection of data representing one state of a roleModel.
#' 
#' @description An S4 class to represent a roleModel state using `localComm`,
#' `metaComm`, & `rolephylo` bojects.
#' 
#' @slot localComm A `localComm` object containing the local community.
#' @slot metaComm A `metaComm` object containing the meta community.
#' @slot phylo A `rolePhylo` object containing the model phylogeny. 
#' 
#' @rdname roleData
#' @include rolePhylo.R
#' @import methods
#' @export

setClass('roleData',
         slots = c(localComm = 'localComm',
                   metaComm = 'metaComm',
                   phylo = 'rolePhylo'))

# constructor
#' @rdname roleData
#' @param localComm A `localComm` object containing the local community.
#' @param metaComm A `metaComm` object containing the meta community.
#' @param phylo A `rolePhylo` object containing the model phylogeny. 
#' @import methods
#' @export

roleData <- function(localComm, metaComm, phylo) {
  return(new('roleData',
             localComm = localComm,
             metaComm = metaComm,
             phylo = phylo))
}

