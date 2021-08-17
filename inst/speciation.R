setGeneric('speciation', function(x, ...) standardGeneric('speciation'))

setMethod('speciation', 'rolePhylo', .specFun)

.specFun <- function(x) {
    n <- x@n

    # add edge
    e <- x@e[, 1:2]
    i <- 3
    j <- which(e[, 2] == i)

    # browser()
    e[e > n] <- e[e > n] + 1

    newNode <- 2 * n + 1
    e <- rbind(e, matrix(c(newNode, newNode,
                           e[j, 2], n + 1),
                         ncol = 2))

    e[j, 2] <- newNode

    # augment edge lengths
    l <- x@e[, 3]
    l <- c(l, 0, 0)
    l[e[, 2] <= n + 1] <- l[e[, 2] <= n + 1] + 1

    # over-write `x` with new info
    x@e <- cbind(e, l)
    x@alive[n + 1] <- TRUE
    x@tipNames[n + 1] <- paste('t', n + 1)
    x@n <- n + 1

    return(x)
}

set.seed(1)
foo <- ape::rphylo(5, 1, 0.1)
foo$edge.length <- foo$edge.length * 5
plot(foo)
axis(1)

boo <- as(foo, 'rolePhylo')
doo <- .specFun(boo)
bla <- as(doo, 'phylo')
bla$edge
