#' @title RoLE phylogeny
#' 
#' @description An S4 class to specify a phylogeny in a way that can be easily
#'     evolved through simulation
#'
#' @slot n number of tips
#' @slot e edge matrix; two columns give ancestor, child pair
#' @slot l numeric vector of edge lengths (in units of time steps = 1/J
#'     generations)
#' @slot alive vector indicating whether tips are extant or not
#' @slot tipNames vector of tip names
#' @slot scale time scale translation to years
#'
#' @export

setClass('rolePhylo',
         slots = c(n = 'numeric',
                   e = 'matrix',
                   l = 'numeric',
                   alive = 'logical',
                   tipNames = 'character',
                   scale = 'numeric'))


#' @title Specify a RoLE model phylogeny
#'
#' @param n number of tips
#' @param e edge matrix; two columns give ancestor, child pair
#' @param l numeric vector of edge lengths (in units of time steps = 1/J
#' generations)
#' @param alive vector indicating whether tips are extant or not
#' @param tipNames vector of tip names
#' @param scale time scale translation to years
#'
#' @export

rolePhylo <- function(n, e, l, alive, tipNames, scale) {
    new('rolePhylo',
        n = n, e = e, l = l, alive = alive, tipNames = tipNames, scale = scale)
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
    
    if(nrow(object@e) != length(object@l)) {
        checks <- c(checks,
                    'unequal number of edges in @e and edge lengths in @l')
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

# register ape phylo
# setOldClass('phylo')

# set coercion method from ape::phylo to roleR::rolePhylo
setAs(from = 'phylo', to = 'rolePhylo',
      def = function(from) {
          # extract number of times
          n <- ape::Ntip(from)
          
          # extract edge matrix and edge lengths
          e <- from$edge
          l <- from$edge.length
          
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
          e <- rbind(e, matrix(-1, nrow = addOn, ncol = 2))
          l <- c(l, rep(0, addOn))
          alive <- c(alive, rep(FALSE, addOn))
          tipNames <- c(tipNames, rep('', addOn))
          
          
          return(rolePhylo(n = n, e = e, l = l, alive = alive,
                           tipNames = tipNames, scale = scale))
      }
)


# set coercion method from roleR::rolePhylo to ape::phylo
setAs(from = 'rolePhylo', to = 'phylo',
      def = function(from) {
          i <- 2 * (from@n - 1)
          
          y <- list(edge = from@e[1:i, ], edge.length = from@l[1:i],
                    tip.label = from@tipNames[1:from@n],
                    Nnode = from@n - 1)
          
          # make any possible 0 or negative edge lengths equal to
          # very small number
          y$edge.length[y$edge.length <= 0] <- .Machine$double.eps
          
          class(y) <- 'phylo'
          
          return(y)
      }
)
