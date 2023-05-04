#' @title A phylogeny of all the species in a `roleData` object.
#' 
#' @description An S4 class to specify a phylogeny for the purpose of the RoLE model.
#'
#' @slot n The number of tips in the phylogeny
#' @slot e The numeric edge matrix of the phylogeny.
#' Each row contains an ancestor-child pair where the 1st column is the ancestor and the 2nd is the child
#' @slot l A numeric vector of edge lengths.
#' The units of l are the time steps (iterations) of the model.
#' Each time step unit is equal to 1/J generations where J is the number of individuals in the local community
#' @slot alive A logical vector indicating whether each tip is extant or not.
#' @slot tipNames A character vector of the tip names.
#' @slot scale A single numeric value of time scale translation to years.
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
#' @param n The number of tips in the phylogeny
#' @param e The numeric edge matrix of the phylogeny.
#' Each row contains an ancestor-child pair where the 1st column is the ancestor and the 2nd is the child
#' @param l A numeric vector of edge lengths.
#' The units of l are the time steps (iterations) of the model.
#' Each time step unit is equal to 1/J generations where J is the number of individuals in the local community
#' @param alive A logical vector indicating whether each tip is extant or not.
#' @param tipNames A character vector of the tip names.
#' @param scale  A single numeric value of time scale translation to years.
#' @return A `rolePhylo` object.
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

