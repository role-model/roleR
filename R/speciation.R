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

.specComm <- function(x) {
    # update number of species
    x@Smax <- x@Smax + 1

    # initialize abundance for new species
    x@abundance[x@Smax] <- 1

    return(x)
}

setMethod('speciation', 'comm', .specComm)


# function to implement speciation for \code{localComm} class objects
#' @param x an object of class \code{localComm}

.specLocal <- function(x, i, p) {

    # update Smax and initialize abundance
    x <- .specComm(x)

    # update abundance
    x@abundance[x@Smax] <- 1

    # index of where unrealized traits begin
    j <- min(which(is.na(x@traits[, 1])))

    # add trait
    x@traits[j, 1] <- x@Smax
    x@traits[j, 2] <- x@traits[i, 2] + rnorm(1, 0, p@params$trait_sigma) # need to figure this out

    return(x)
}
setMethod('speciation', 'localComm', .specLocal)


#' function to implement speciation for \code{roleModel} class objects
#' @param x an object of class \code{roleModel}
#'
.specRole <- function(x) {

    # note: `Smax` from `@localComm` and `@metaComm` and `n` from `@phylo` are
    # all enforced to be equal, so we can sample from any but we have to
    # weight the probabilities by abundances and immigration

    # dispersal prob
    dp <- x@params@params$dispersal_prob

    # normalized abundances at meta and local levels
    mp <- x@metaComm@abundance[1:x@metaComm@Smax]
    mp <- mp / sum(mp)
    lp <- x@localComm@abundance[1:x@localComm@Smax]
    lp <- lp / sum(lp)

    # prob of selecting a parent
    # metacomm abundance weighted by dispersal prob + local comm abundance weighted by birth
    pp <- dp * mp + (1 - dp) * lp

    # index of parent
    i <- sample(x@phylo@n, size = 1, prob = pp)

    # update slots of the role model object
    x@localComm <- speciation(x@localComm, i = i, p = x@params)
    x@metaComm <- speciation(x@metaComm)
    x@phylo <- speciation(x@phylo, i = i)

    return(x)
}
setMethod('speciation', 'roleModel', .specRole)
