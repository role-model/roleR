#' @title Death
#'
#' @description Generic and methods for computing the death process
#'
#' @param x the object which determines method dispatch
#' @param i the index of the species undergoing death
#'
#' @rdname death
#' @export

setGeneric('death',
           function(x, i, ...) standardGeneric('death'),
           signature = 'x')

#' function to implement death for \code{*comm} class objects
#' @param x an object of class \code{localComm}
#' @param i the index of the species undergoing death
#' @param params a \code{roleParams} object

.deathComm <- function(x, i) {

    # decrement abundance
    x@abundance[i] <- x@abundance[i] - 1

    # if abundance is zero, make species extinct (?)
    if(x@abundance[i] <= 0)
    {
        #does comm need to hold a vector of species that are extinct?
        #to prevent births, immigrations, etc for extinct
    }

    return(x)
}

setMethod('death', 'comm', .deathComm)

#' function to implement death for \code{*rolePhylo} class objects
#' @param x an object of class \code{rolePhylo}
#' @param i the index of the tip undergoing death
#' @param params a \code{roleParams} object

.deathPhylo <- function(x, i) {

    # find tips where species abundance is <=0 in iterSim and pass as i
    # or maybe do this as part of .deathPhylo?
    x@alive[i] <- FALSE
    # do tip indices follow species indices?

    return(x)
}

setMethod('death', 'rolePhylo', .deathPhylo)
