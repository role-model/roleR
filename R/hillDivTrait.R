# `hillDiv` method for class `trait`

hillDivTrait <- function(x, traits, q) {
    p <- x / sum(x)
    dij <- as.matrix(dist(traits))
    Q <- as.vector(p %*% dij %*% p)
    a <- outer(p, p, '*') / Q

    Hk <- sapply(q, function(qk) {
        if(qk == 1) {
            return(exp(- sum(dij * a * log(a))))
        } else {
            return(sum(dij * a^qk)^(1 / (1 - qk)))
        }
    })

    D <- sqrt(Hk / Q)

    return(D)
}



# test (needs dev version of hillR)
# fooComm <- matrix(1:10, nrow = 1)
# colnames(fooComm) <- letters[1:10]
# fooTraits <- data.frame(traits = 1:10, row.names = letters[1:10])
# fooTraits[colnames(fooComm), , drop = FALSE]
#
#
# hillR::hill_func(fooComm, fooTraits, q = 1.5)
#
# hillDivTrait(1:10, 1:10, 1.5)
