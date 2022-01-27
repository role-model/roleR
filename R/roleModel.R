# This contains R wrappers for roleModel and containing C++ objects

setClass('localComm',
         slots = c(abundanceIndv = 'localComm',
                   metaComm = 'metaComm',
                   phylo = 'rolePhylo',
                   params = 'roleParams'))

localComm <- function(local, meta, phy, p) {
  new('localComm',
      localComm = local,
      metaComm = meta,
      phylo = phy,
      params = p)s
}


#' @title An S4 class to specify the entire RoLE model
#'
#' @slot localComm an object of class \code{localComm}
#' @slot metaComm an object of class \code{metaComm}
#' @slot phylo an object of class \code{rolePhylo}
#' @slot params an object of calss \code{roleParams}
#'
#' @export
#' @include comm.R roleParams.R rolePhylo.R


setClass('roleModel',
         slots = c(localComm = 'localComm',
                   metaComm = 'metaComm',
                   phylo = 'rolePhylo',
                   params = 'roleParams'))


#' @title Create a RoLE model object
#'
#' @param localComm an object of class \code{localComm}
#' @param metaComm an object of class \code{metaComm}
#' @param phylo an object of class \code{rolePhylo}
#' @param params an object of class \code{roleParams}
#'
#' @return an object of class \code{roleModel}
#'
#' @export


roleModel <- function(local, meta, phy, p) {
  new('roleModel',
      localComm = local,
      metaComm = meta,
      phylo = phy,
      params = p)
}
