#' @title A phylogeny of all the species of a `roleData`
#' 
#' @description An S4 class to specify a phylogeny for the purpose of the RoLE model
#'
#' @slot n the number of tips in the phylogeny
#' @slot e the numeric edge matrix of the phylogeny.
#' Each row contains an ancestor-child pair where the 1st column is the ancestor and the 2nd is the child
#' @slot l a numeric vector of edge lengths.
#' The units of l are the time steps (iterations) of the model.
#' Each time step unit is equal to 1/J generations where J is the number of individuals in the local community
#' @slot alive a logical vector indicating whether each tips is extant or not
#' @slot tipNames a character vector of the names of each tip 
#' @slot scale a single numeric value of time scale translation to years
#' 
#' @rdname rolePhylo
#' @export

setClass('rolePhylo',
         slots = c(n = 'numeric',
                   e = 'matrix',
                   l = 'numeric',
                   alive = 'logical',
                   tipNames = 'character',
                   scale = 'numeric'))

#' @title Create a rolePhylo
#'
#' @param n 
#' @param e 
#' @param l 
#' @param alive 
#' @param tipNames 
#' @param scale
#' @return a `rolePhylo` object
#' 
#' @rdname rolePhylo
#' @export

rolePhylo <- function(n, e, l, alive, tipNames, scale) {
    new('rolePhylo',
        n = n, e = e, l = l, alive = alive, tipNames = tipNames, scale = scale)
}

# checker function for validation
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
