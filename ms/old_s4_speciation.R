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