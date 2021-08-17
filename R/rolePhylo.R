#' @title An S4 class to specify a phylogeny in a way that can be easily
#' evolved through simulation
#'
#' @slot n number of tips
#' @slot e extended edge matrix; first two columns give ancestor, child pair,
#' their column is edge lengths (in units of time steps = 1/J generations)
#' @slot alive vector indicating whether tips are extant or not
#' @slot tipNames vector of tip names
#' @slot scale time scale translation to years
#'
#' @export

setClass('rolePhylo',
         slots = c(n = 'numeric',
                   e = 'matrix',
                   alive = 'logical',
                   tipNames = 'character',
                   scale = 'numeric'))


#' @title Specify a RoLE model phylogeny
#'
#' @param n number of tips
#' @param e extended edge matrix; first two columns give ancestor, child pair,
#' their column is edge lengths (in units of time steps = 1/J generations)
#' @param alive vector indicating whether tips are extant or not
#' @param tipNames vector of tip names
#' @param scale time scale translation to years
#'
#' @export

rolePhylo <- function(n, e, alive, tipNames, scale) {
    new('rolePhylo',
        n = n, e = e, alive = alive, tipNames = tipNames, scale = scale)
}

# checker function for validation
#' @param object an object of class localComm

checkRolePhylo <- function(object) {
    checks <- c()

    if(length(object@alive) < length(object@tipNames)) {
        checks <- c(checks, 'not all named tips are represented in @alive')
    }

    if(length(object@alive) < object@n) {
        checks <- c(checks,
                    'fewer tips specified in @alive than indicated by @n')
    }

    if(length(object@tipNames) < object@n) {
        checks <- c(checks,
                    'fewer tips specified in @tipNames than indicated by @n')
    }

    if(nrow(object@e) < 2 * (object@n - 1)) {
        checks <- c(checks,
                    'edge matrix does not contain sufficient rows for number of tips')
    }

    # if any issues, return them, otherwise all OK
    if(length(checks) > 0) {
        return(checks)
    } else {
        return(TRUE)
    }
}


# validate
setValidity('rolePhylo', checkRolePhylo)


# set coercion method from ape::phylo to roleR::rolePhylo
setAs(from = 'phylo', to = 'rolePhylo',
      def = function(from) {
          # extract number of times
          n <- ape::Ntip(from)

          # extract edge matrix and edge lengths
          e <-  cbind(from$edge, from$edge.length)

          # extract tip labels
          tipNames <- from$tip.label

          # calculate alive or not

          tipAge <- ape::node.depth.edgelength(from)[1:n]

          alive <- rep(TRUE, n)
          alive[tipAge < max(tipAge)] <- FALSE


          # set default scale
          scale <- 1


          # buffer objects so we can add new species without augmenting objects
          addOn <- n * 100
          e <- rbind(e, matrix(NA, nrow = addOn, ncol = 3))
          alive <- c(alive, rep(FALSE, addOn))


          return(rolePhylo(n = n, e = e, alive = alive,
                           tipNames = tipNames, scale = scale))
      }
)


# set coercion method from roleR::rolePhylo to ape::phylo
setAs(from = 'rolePhylo', to = 'phylo',
      def = function(from) {
          i <- 2 * (from@n - 1)

          y <- list(edge = from@e[1:i, 1:2], edge.length = from@e[1:i, 3],
                    tip.label = from@tipNames,
                    Nnode = from@n - 1)

          class(y) <- 'phylo'

          return(y)
      }
)

# test
# tre <- ape::rphylo(5, 1, 0.1)
# plot(tre)
# ape::nodelabels()
# foo <- as(tre, 'rolePhylo')
# boo <- as(foo, 'phylo')
# plot(boo)
