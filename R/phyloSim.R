#' @title phylo simulation
#'
#' @param x tree
#' @param i index
#'
#' @return phylo
#'
#' @export

phyloSim <- function(x, i) {
    newX <- as(speciation(as(x, 'rolePhylo'), i), 'phylo')

    return(newX)
}

# test
# set.seed(1)
# tre <- ape::rphylo(5, 1, 0.1)
# plot(tre)
# axis(1)
#
# newTre <- phyloSim(tre, 3)
# split.screen(c(1, 1), erase = FALSE)
# plot(newTre, edge.color = 'red')
# axis(1, col = 'red')
#
# speciation(as(tre, 'rolePhylo'), 3)
