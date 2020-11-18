#' @title Hill number diversities
#'
#' @description Hill number diversities for species abundance, trait, population
#' genetic, and phylogenetic data
#'
#'
#' @param x a vector of species abundances
#' @param q a vector of Renyei entropy scales
#'
#' @details Stub
#'
#' @examples
#' hillDivSAD(1:10, 0:3)
#' hillDivTrait(1:10, 1:10, 0:3)
#'
#' @return A vector, same length as \code{q}, giving diversities for each
#' value of \code{q}
#'
#' @rdname hillDiv
#' @export

hillDivSAD <- function(x, q) {
    x <- x/sum(x)

    hill <- outer(x, q, '^')
    hill <- colSums(hill)^(1 / (1 - q))
    hill[hill == 1] <- exp(sum(-x * log(x)))

    return(hill)
}


# test
# hillDivSAD(1:10, 0:3)
# sapply(0:3, function(qi) hillR::hill_taxa(1:10, qi))




#' @rdname hillDiv
#' @export
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
