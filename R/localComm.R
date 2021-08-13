#' @title An S4 class to specify the local community of the RoLE model
#'
#' @slot abundance a numeric vector of abundances for each species
#' @slot traits matrix of traits; the first column specifies the species index
#' (i.e. the index of that species in the \code{abundance} vector) and the
#' second column specifies the trait value
#' @slot pi a numeric vector of genetic diversities for each species
#' @slot Smax a single integer specifying the total number of species ever
#' recorded in the local community (both locally extinct and extant)
#'
#' @export

setClass('localComm',
         slots = c(abundance = 'numeric',
                   traits = 'matrix',
                   pi = 'numeric',
                   Smax = 'numeric'))


#' @title Specify a local community object
#'
#' @param abundance a numeric vector of abundances for each species
#' @param traits matrix of traits; the first column specifies the species index
#' (i.e. the index of that species in the \code{abundance} vector) and the
#' second column specifies the trait value
#' @param pi a numeric vector of genetic diversities for each species
#' @param Smax a single integer specifying the total number of species ever
#' recorded in the local community (both locally extinct and extant)
#'
#' @details \code{Smax} is used for internal bookkeeping.  The dimensions of
#' \code{abundance}, \code{traits}, and \code{pi} can be greater than
#' \code{Smax}.  In such cases, \code{Smax} is used to index where to begin
#' adding new information (e.g. a new species can be added at the index
#' \code{Smax + 1}).
#'
#' @return an object of class \code{localComm}
#'
#' @seealso \code{\link{localComm-class}}
#'
#' @export

localComm <- function(abundance, traits, pi) {
    new('localComm',
        abundance = abundance, traits = traits, pi = pi)
}


#' checker function for validation
#' @param object an object of class localComm

checkLocalComm <- function(object) {
    checks <- c()

    if(length(object@abundance) != length(object@pi)) {
        checks <- c(checks,
                    'lenths of @abundance and @pi must be equal')
    } else {
        if(length(object@abundance) < object@Smax) {
            checks <- c(checks,
                        '@Smax greater than number of species with data')
        }
    }

    if(max(object@traits[, 1]) > object@Smax) {
        checks <- c(checks,
                    '@traits records traits for more species than allowed by @Smax')
    }

    # if any issues, return them, otherwise all OK
    if(length(checks) > 0) {
        return(checks)
    } else {
        return(TRUE)
    }
}


#' validate
setValidity('localComm', checkLocalComm)
