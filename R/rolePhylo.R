#' @title An S4 class to specify a phylogeny 
#'
#' @param ntips number of tips
#' @param edges edge matrix; two columns give ancestor, child pair
#' @param lengths numeric vector of edge lengths (in units of time steps = 1/J
#' generations)
#' @param alive vector indicating whether tips are extant or not
#' @param tipNames vector of tip names
#' @param scale time scale translation to years
#'
#' @export

setClass('rolePhylo',
         slots = c(ntips = 'numeric',
                   edges = 'matrix',
                   lengths = 'numeric',
                   alive = 'logical',
                   tipNames = 'character',
                   scale = 'numeric'), 
         validity=function(object)
         {
           checks <- c()
           
           if(length(object@alive) < length(object@tipNames)) {
             checks <- c(checks, 'not all named tips are represented in @alive')
           }
           
           if(length(object@alive) < object@ntips) {
             checks <- c(checks,
                         'fewer tips specified in @alive than indicated by @n')
           }
           
           if(length(object@tipNames) < object@ntips) {
             checks <- c(checks,
                         'fewer tips specified in @tipNames than indicated by @n')
           }
           
           if(nrow(object@edges) < 2 * (object@ntips - 1)) {
             checks <- c(checks,
                         'edge matrix does not contain sufficient rows for number of tips')
           }
           
           if(nrow(object@edges) != length(object@lengths)) {
             checks <- c(checks,
                         'unequal number of edges in @e and edge lengths in @l')
           }
           
           # if any issues, return them, otherwise all OK
           if(length(checks) > 0) {
             return(checks)
           } else {
             return(TRUE)
           }
         })


#' @title Specify a RoLE model phylogeny
#'
#' @param ntips number of tips
#' @param edges edge matrix; two columns give ancestor, child pair
#' @param lengths numeric vector of edge lengths (in units of time steps = 1/J
#' generations)
#' @param alive vector indicating whether tips are extant or not
#' @param tipNames vector of tip names
#' @param scale time scale translation to years
#'
#' @export

rolePhylo <- function(ntips, edges, lengths, alive, tipNames, scale) {
  new('rolePhylo',
      ntips = ntips, edges = edges, lengths = lengths, alive = alive, tipNames = tipNames, scale = scale)
}

rolePhyloToCpp <- function(phylo){
  n <- phylo@ntips
  e <- phylo@edges
  l <- phylo@lengths
  alive <- phylo@alive
  tipNames <- phylo@tipNames
  scale <- phylo@scale
  out <- new(rolePhyloCpp,n,e,l,alive,tipNames,scale)
  return(out)
}

rolePhyloFromCpp <- function(phylo){
  n <- phylo$n
  e <- phylo$e
  l <- phylo$l
  alive <- phylo$alive
  tipNames <- phylo$tipNames
  scale <- phylo$scale
  out <- rolePhylo(n,e,l,alive,tipNames,scale)
  return(out)
}