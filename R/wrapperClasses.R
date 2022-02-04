# wrapperClasses.R contains simple S4 object wrappers for roleModel constituents

#' @title An S4 class to specify a local community within a ROLE model 
#'
#' @slot abundanceIndv a numeric vector of 0s and 1s specifying the binary abundance (alive or dead) 
#' of each individual
#' @slot speciesIDsIndv a numeric vector specifying species of each individual as
#' an index of the species vectors
#' @slot traitsIndv a numeric vector of trait values for each individual
#' @slot abundanceSp abundances of every species i held at index i
#' @slot traitsSp traits of every species i held at index i
#' @slot sequencesSp genetic diversities of species
#' @export

setClass('localComm',
         slots = c(abundanceIndv = 'numeric',
                   speciesIDsIndv = 'numeric',
                   traitsIndv = 'numeric',
                   abundanceSp = 'numeric',
                   traitsSp = 'numeric',
                   sequencesSp = 'numeric'))
# constructor 
localComm <- function(abundanceIndv, speciesIDsIndv, traitsIndv, 
                      abundanceSp, traitsSp,sequencesSp) {
  return(new('localComm',
      abundanceIndv = abundanceIndv,
      speciesIDsIndv = speciesIDsIndv,
      traitsIndv = traitsIndv,
      abundanceSp = abundanceSp,
      traitsSp = traitsSp,
      sequencesSp = sequencesSp))
}

#' @title An S4 class to specify a meta community within a ROLE model 
#'
#' @slot abundanceSp abundances of every species i held at index i
#' @slot traitsSp traits of every species i held at index i
#' @export

setClass('metaComm',
         slots = c(abundanceSp = 'numeric',
                   traitsSp = 'numeric'))
# constructor
metaComm <- function(abundanceSp,traitsSp,Smax) {
  return(new('metaComm',
      abundanceSp = abundanceSp,
      traitsSp = traitsSp))
}

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
                   scale = 'numeric'))


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
      ntips = ntips, e = e, l = l, alive = alive, tipNames = tipNames, scale = scale)
}

#' @title An S4 class to role model data for timeseries
#' @slot localComm an object of class \code{localComm}
#' @slot metaComm an object of class \code{metaComm}
#' @slot phylo an object of class \code{rolePhylo}
#' @slot stats an object of class \code{data-frame} - 1st col contains param names,
#' 2nd col contains numeric entropies used to compute hill statistics, 
#' and 3rd col contains hill statistic values
#' 
#' @export

setClass('roleData',
         slots = c(localComm = 'localComm',
                   metaComm = 'metaComm',
                   phylo = 'rolePhylo',
                   stats = 'data.frame'))

# constructor
roleData <- function(localComm,metaComm,phylo) {
  return(new('roleData',
      localComm = localComm,
      metaComm = metaComm,
      phylo = phylo))
}