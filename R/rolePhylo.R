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
         slots = c(abundance = 'numeric',
                   traits = 'matrix',
                   pi = 'numeric',
                   Smax = 'numeric'))


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

rolePhylo <- function() {
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
