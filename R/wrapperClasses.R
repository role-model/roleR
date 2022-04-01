# wrapperClasses.R contains ple S4 object wrappers for roleModel constituents

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
         slots = c(abundanceSp = 'numeric',
                   traitsSp = 'matrix'))
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
                   stats = 'data.frame',
                   iterNum = 'integer'))

# constructor
roleData <- function(localComm,metaComm,phylo,iterNum) {
  return(new('roleData',
      localComm = localComm,
      metaComm = metaComm,
      phylo = phylo,
      iterNum = iterNum))
}

roleDataFromCpp <- function(data) {
  local <- localComm(data$local$abundance_indv,data$local$species_ids,
                     data$local$traits,data$local$abundance_sp,
                     data$local$traits_sp,data$local$pi_sp)
  
  meta <- metaComm(data$meta$abundance,data$meta$traits)
  phylo <- rolePhyloFromCpp(data$phylo)
  stats <- data$stats 
  iterNum <- data$iter_num
  
  return(roleData(local,meta,phylo,iterNum))
}