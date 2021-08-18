# private generic for speciation
#' @param x the object which determines method dispatch
#' @param i the index of the tip undergoing speciation

setGeneric('speciation',
           function(x, i) standardGeneric('speciation'),
           signature = 'x')


# function to implement speciation for \code{rolePhylo} class objects
#' @param x an object of class \code{rolePhylo}
#' @param i the index of the tip undergoing speciation
#' @note method is set below after the function def

.specFun <- function(x, i) {
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


# set the (private) method
setMethod('speciation', 'rolePhylo', .specFun)

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
