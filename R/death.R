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

#' function to implement death for \code{*roleModel} class objects
#' @param x an object of class \code{roleModel}
#' @param i the index of the species undergoing death

.deathRoleModel <- function(x) {

    # sample a species for death proportional to species abundance
    i <- sample(x@localComm@Smax, 1, prob = x@localComm@abundance[1:x@localComm@Smax])

    # call death on localComm
    x@localComm <- death(x@localComm, i)

    # if death led to extinction, call death on rolePhylo
    if(x@localComm@abundance[i] <= 0)
        x@rolePhylo <- death(x@rolePhylo, i)

    return(x)
}

setMethod('death', 'roleModel', .deathRoleModel)

#' function to implement death for \code{*comm} class objects
#' @param x an object of class \code{localComm}
#' @param i the index of the species undergoing death

.deathComm <- function(x, i) {

    # decrement abundance
    x@abundance[i] <- x@abundance[i] - 1

    return(x)
}

setMethod('death', 'comm', .deathComm)

#' function to implement death for \code{*rolePhylo} class objects
#' @param x an object of class \code{rolePhylo}
#' @param i the index of the tip undergoing death
#' @param params a \code{roleParams} object

.deathPhylo <- function(x, i) {

    # set tip to dead
    x@alive[i] <- FALSE

    return(x)
}

setMethod('death', 'rolePhylo', .deathPhylo)

# TEST
# source("R/comm.R")
# library("ape")
# source("R/rolePhylo.R")
# source("R/roleParams.R")
# source("R/roleModel.R")
# source("R/death.R")
# foo <- localComm(1:10,matrix(1:100,nrow = 10, ncol = 10),1:10,10)
# doo <- ape::rphylo(5, 1, 0.1)
# doo <- as(doo, 'rolePhylo')
# #todo finish phylo tests
# foo@abundance[1]
# foo <- death(foo,1)
# foo@abundance[1]

