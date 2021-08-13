#' @title An S4 class to specify the entire hierrarchical community of the
#' RoLE model
#'
#' @slot localComm an object of class \code{localComm}
#' @slot params an object of calss \code{roleParams}
#'
#' @export
#' @include localComm.R roleParams.R


setClass('roleComm',
         slots = c(localComm = 'localComm',
                   # metaComm = 'metaComm',
                   # phylo = 'rolePhylo',
                   params = 'roleParams'))
