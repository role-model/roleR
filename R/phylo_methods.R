#' make a class `rolePhylo` that has
#' @slot n number of tips
#' @slot e extended edge matrix; first two columns give ancestor, child pair,
#' their column is edge lengths (in units of time steps = 1/J generations)
#' @slot alive vector indicating whether tips are extant or not
#' @slot tipNames vector of tip names
#' @slot scale time scale translation to years


# replace this eventually with `as`, see example here:
# https://github.com/acidgenomics/r-pointillism/blob/main/R/coerce-methods.R

as_rolePhylo <- function(x, alive = NULL, scale = NULL) {
    # extract number of times
    n <- ape::Ntip(x)

    # extract edge matrix and edge lengths
    e <-  cbind(x$edge, x$edge.length)

    # extract tip labels
    tipNames <- x$tip.label

    # calculate alive or not if necessary
    if(is.null(alive)) {
        tipAge <- ape::node.depth.edgelength(x)[1:n]

        alive <- rep(TRUE, n)
        alive[tipAge < max(tipAge)] <- FALSE
    }

    # set default scale if necessary
    if(is.null(scale)) {
        scale <- 1
    }

    # buffer objects so we can add new species without augmenting objects
    addOn <- n * 100
    e <- rbind(e, matrix(NA, nrow = addOn, ncol = 2))
    alive <- c(alive, rep(FALSE, addOn))

    # output
    out <- list(n = n,
                e = e,
                alive = alive,
                tipNames = tipNames,
                scale = scale)

    class(out) <- 'rolePhylo'

    return(out)
}
