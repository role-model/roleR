#' @title Speciation
#'
#' @description Generic and methods for computing the speciation process
#'
#' @param x the object which determines method dispatch
#' @param i the index of the species undergoing speciation
#'
#' @rdname speciation
#' @export

setGeneric('speciation',
           function(x, i, ...) standardGeneric('speciation'),
           signature = 'x')


# function to implement speciation for \code{rolePhylo} class objects
#' @param x an object of class \code{rolePhylo}
#' @param i the index of the tip undergoing speciation
#' @param ... additional parameters passed to specific methods
#' @note method is set below after the function def

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

# test
# foo <- ape::rphylo(5, 1, 0.1)
# foo$edge.length <- foo$edge.length * 5
# plot(foo)
# axis(1)
#
# boo <- as(foo, 'rolePhylo')
#
# doo <- speciation(boo, 3)
# bla <- as(doo, 'phylo')
# plot(bla)


# function to implement speciation for \code{*comm} class objects
#' @param x an object of class \code{localComm}
#' @param i the index of the species undergoing speciation

.specComm <- function(x, i) {
    # update number of species
    x@Smax <- x@Smax + 1

    # add abundance to new species
    x@abundance[x@Smax] <- 1

    return(x)
}

setMethod('speciation', 'comm', .specComm)


# function to implement speciation for \code{localComm} class objects
#' @param params a \code{roleParams} object

.specLocal <- function(x, i, params) {
    # update abund and Smax
    x <- .specComm(x, i)

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

# test
# foo <- c(1:5, rep(0, 5))
# x <-  localComm(foo, cbind(foo, foo), foo, 5)
# p <- roleParams(list(trait_sigma = 0.0001), 'sim')
# speciation(x, 5, params = p)
