#' @title An S4 class to specify a local community within a ROLE model 
#'
#' @slot abundanceIndv a numeric vector of 0s and 1s specifying the binary abundance (alive or dead) 
#' of each individual
#' @slot speciesIDsIndv a numeric vector specifying species of each individual as
#' an index of the species vectors
#' @slot traitsIndv a numeric vector of trait values for each individual
#' @slot abundanceSp abundances of every species i held at index i
#' @slot traitsSp traits of every species i held at index i
#' @slot gdiversitiesSp genetic diversities of species
#' @export

setClass('localComm',
         slots = c(abundanceIndv = 'numeric',
                   speciesIDsIndv = 'numeric',
                   traitsIndv = 'numeric',
                   abundanceSp = 'numeric',
                   traitsSp = 'numeric',
                   gdiversitiesSp = 'numeric'))
# constructor 
localComm <- function(abundanceIndv, speciesIDsIndv, traitsIndv, 
                      abundanceSp, traitsSp,gdiversitiesSp) {
  return(new('localComm',
      abundanceIndv = abundanceIndv,
      speciesIDsIndv = speciesIDsIndv,
      traitsIndv = traitsIndv,
      abundanceSp = abundanceSp,
      traitsSp = traitsSp,
      gdiversitiesSp = gdiversitiesSp))
}

#' @title An S4 class to specify a meta community within a ROLE model 
#'
#' @slot abundanceSp abundances of every species i held at index i
#' @slot traitsSp traits of every species i held at index i
#' @export

setClass('metaComm',
         slots = c(abundanceSpM = 'numeric',
                   traitsSpM = 'matrix'))
# constructor
metaComm <- function(abundanceSpM,traitsSpM) {
  return(new('metaComm',
      abundanceSpM = abundanceSpM,
      traitsSpM = traitsSpM))
}