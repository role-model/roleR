#' @title Speciation
#'
#' @description Generic and methods for computing the speciation process
#'
#' @param x the object which determines method dispatch
#' @param ... additional parameters passed to specific methods
#'
#' @rdname speciation
#' @export

setGeneric('speciation',
           function(x, ...) standardGeneric('speciation'),
           signature = 'x')


# function to implement speciation for \code{rolePhylo} class objects
#' @param x an object of class \code{rolePhylo}
#' @param i the index of the tip undergoing speciation

.specPhylo <- function(x, i) {
    # number of tips
    n <- x@n

    # index of where unrealized edges in edge matrix start
    eNew <- min(which(x@e[, 1] == -1))

    # index of where to add new edge
    j <- which(x@e[, 2] == i)


    # add one to internal node indices
    x@e[x@e > n] <- x@e[x@e > n] + 1

    # add new node
    newNode <- 2 * n + 1 # index of new node
    x@e[(0:1) + eNew, 1] <- newNode # add internal node
    x@e[eNew, 2] <- x@e[j, 2] # add old tip
    x@e[eNew + 1, 2] <- n + 1 # add new tip

    # add edge connecting parent to new node
    x@e[j, 2] <- newNode

    # augment edge lengths
    x@l[(0:1) + eNew] <- 0 # add new edges
    # x@l[x@e[, 2] <= n + 1] <- x@l[x@e[, 2] <= n + 1] + 1 # increase tip length

    # over-write other slots of `x` with new info
    x@alive[n + 1] <- TRUE
    x@tipNames[n + 1] <- paste0('t', n + 1)
    x@n <- n + 1

    return(x)
}


# set the method
setMethod('speciation', 'rolePhylo', .specPhylo)

# TEST
# source("R/speciation.R")
# foo <- ape::rphylo(5, 1, 0.1)
# foo$edge.length <- foo$edge.length * 5
# plot(foo)
# axis(1)
#
# boo <- as(foo, 'rolePhylo')
#
# doo <- speciation(boo, i = 3)
# bla <- as(doo, 'phylo')
# plot(bla)


# function to implement speciation for \code{*comm} class objects
#' @param x an object of class \code{comm}
#' @param i the index of the species undergoing speciation

.specComm <- function(x, i) {
    # update number of species
    x@Smax <- x@Smax + 1

    # initialize abundance for new species
    x@abundance[x@Smax] <- 0

    return(x)
}

setMethod('speciation', 'comm', .specComm)


# function to implement speciation for \code{localComm} class objects
#' @param x an object of class \code{localComm}
#' @param i the index of the species undergoing speciation
#' @param params a \code{roleParams} object

.specLocal <- function(x, i, params) {
    # update Smax and initialize abundance
    x <- .specComm(x, i)

    # update abundance
    x@abundance[x@Smax] <- 1

    # index of where unrealized traits begin
    # note: we need to do `Smax - 1` because above where we did
    # `.specComm(x, i)` that already updated `Smax`
    j <- max(which(x@traits[, 1] == x@Smax - 1)) + 1

    # add trait
    x@traits[j, 1] <- x@Smax
    x@traits[j, 2] <- x@traits[i, 2] + rnorm(1, 0, params@params$trait_sigma) # need to figure this out

    return(x)
}

setMethod('speciation', 'localComm', .specLocal)

#' function to implement speciation for \code{roleModel} class objects
#' @param params a \code{roleParams} object

.specRoleModel <- function(x, params) {

    # add sampling for i
    x@localComm <- speciation(x@localComm)
    x@rolePhylo <- speciation(x@rolePhylo)

    return(x)
}

setMethod('speciation', 'localComm', .specLocal)