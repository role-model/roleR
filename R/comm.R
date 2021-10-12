#' @title S4 classes to specify community objects
#'
#' @description The \code{comm} class is used as the parent class for
#' \code{localComm} and \code{metaComm} child classes.
#'
#' @slot abundance a numeric vector of abundances for each species
#' @slot traits matrix of traits; the first column specifies the species index
#' and the subsequent columns specify the trait values
#' @slot Smax a single integer specifying the total number of species ever
#' recorded in the local community (both locally extinct and extant)
#'
#' @rdname comm-class
#' @export

setClass('comm',
         slots = c(abundance = 'numeric',
                   traits = 'matrix',
                   Smax = 'numeric'))


# Local Community
#' @slot pi a numeric vector of genetic diversities for each species
#'
#' @rdname comm-class
#' @export

setClass('localComm',
         contains = 'comm',
         slots = c(pi = 'numeric'))

# Meta Community
#'
#' @rdname comm-class
#' @export

setClass('metaComm',
         contains = 'comm')


#' @title Specify community objects
#'
#' @param abundance a numeric vector of abundances for each species
#' @param traits matrix of traits; the first column specifies the species index
#' (i.e. the index of that species in the \code{abundance} vector) and the
#' subsequent columns specify the trait values
#' @param Smax a single integer specifying the total number of species ever
#' recorded in the local community (both locally extinct and extant)
#'
#' @details \code{Smax} is used for internal bookkeeping.  The dimensions of
#' \code{abundance} and \code{traits} can be greater than
#' \code{Smax}.  In such cases, \code{Smax} is used to index where to begin
#' adding new information (e.g. a new species can be added at the index
#' \code{Smax + 1}).
#'
#' @return an object of class \code{comm}
#'
#' @seealso \code{\link{comm-class}}
#'
#' @rdname comm
#' @export

comm <- function(abundance, traits, Smax) {
    new('comm',
        abundance = abundance, traits = traits, Smax = Smax)
}

# Local Community
#' @param pi a numeric vector of genetic diversities for each species
#'
#' @rdname comm
#' @export

localComm <- function(abundance, traits, pi, Smax) {
    new('localComm',
        abundance = abundance, traits = traits, pi = pi, Smax = Smax)
}

# Meta Community
#'
#' @rdname comm
#' @export

metaComm <- function(abundance, traits, Smax) {
    new('metaComm',
        abundance = abundance, traits = traits, Smax = Smax)
}


# checker function for parent comm class validation
#' @param object an object of class comm

checkComm <- function(object) {

    checks <- c()
    if(length(object@abundance) < object@Smax) {
        checks <- c(checks,
                    '@Smax greater than number of species with data')
    }

    if(max(object@traits[, 1]) > object@Smax) {
        checks <- c(checks,
                    '@traits records traits for more species than allowed by @Smax')
    }

    if(nrow(object@traits) != length(object@abundance)) {
        checks <- c(checks,
                    'unequal lengths/nrows of @abundance and @traits')
    }

    # if any issues, return them, otherwise all OK
    if(length(checks) > 0) {
        return(checks)
    } else {
        return(TRUE)
    }
}


# validate parent comm class
#setValidity('comm', checkComm)


# checker function for localComm class validation
#' @param object an object of class localComm

checkLocalComm <- function(object) {
    checks <- c()

    if(length(object@abundance) != length(object@pi)) {
        checks <- c(checks,
                    'lengths of @abundance and @pi must be equal')
    }

    #make sure all indices of trait values match indices of nonzero abundance
    indices = na.omit(object@traits[,1])
    #if all indices with trait values have <= 0 abundance
    if(all(object@abundance[indices] <= 0)){
        checks <- c(checks, 'all indices of trait values match indices of
                    nonzero abundance')
    }

    #make sure local indices dont exceed Smax
    traits = na.omit(object@traits)
    # if traits is not null
    # if last index of traits > Smax
    if(isTRUE(traits)){
        if(tail(traits[,1], n = 1) > object@Smax){
            checks <- c(checks, 'local indices cannot exceed Smax')
        }
    }

    # if any issues, return them, otherwise all OK
    if(length(checks) > 0) {
        return(checks)
    } else {
        return(TRUE)
    }
}


# validate localComm class
#setValidity('localComm', checkLocalComm)
