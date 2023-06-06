#' @title The local community of a `roleData` 
#' 
#' @description An S4 class to specify the state of a local community
#' @slot indSpecies a numeric vector of the species IDs for each individual
#' @slot indTrait a numeric vector of the trait values for each individual
#' @slot indSeqs a character vector of the gene sequences for each individual
#' @slot spGenDiv a numeric vector of the genetic diversities for each species
#' @slot spTrait a numeric vector of the mean trait value of each species
#' @slot spAbund a numeric vector of the abundance (number of individuals) of 
#'     each species
#' @slot spAbundHarmMean numeric vector of the harmonic mean of species 
#'     abundances
#' @slot spLastOriginStep numeric vector holding time of most recent origin of 
#'     each species in the local community
#' @slot spExtinctionStep numeric vector of most recent extirpation step of each 
#'     species in the local community
#' @slot equilibProp numeric proportion of equilibrium achived 
#' @rdname localComm
#' @import methods
#' @export

setClass('localComm',
         slots = c(indSpecies = 'numeric',
                   indTrait = 'numeric',
                   indSeqs = 'character',
                   spGenDiv = 'numeric',
                   spTrait = 'numeric',
                   spAbund = 'numeric',
                   spAbundHarmMean = 'numeric',
                   spLastOriginStep = 'numeric',
                   spExtinctionStep = 'numeric',
                   equilibProp = 'numeric'
         ))


# constructor 
#' @rdname localComm
#' @param indSpecies indSpecies
#' @param indTrait indTrait
#' @param indSeqs indSeqs
#' @param spGenDiv spGenDiv
#' @import methods
#' @export

localComm <- function(indSpecies, indTrait, indSeqs) { 
    
    # create vectors indexed by spp ID from individual-indexed vectors
    spAbund <- as.vector(tabulate(indSpecies))
    spTrait <- as.vector(tapply(indTrait, indSpecies, mean))
    
    # initialize harmonic mean, with only 1 timestep harm mean = abund
    spAbundHarmMean <- spAbund
    
    # origin index is all -1 (i.e. not there yet)
    spLastOriginStep <- rep(-1, length(spAbundHarmMean))
    
    # except change origin to 0 for those that are there
    spLastOriginStep[spAbund > 0] <- 0 
    
    # no one extinct yet so all -1
    spExtinctionStep <- rep(-1, length(spAbundHarmMean))
    
    # all 0 for gen div (will be over-written by msprime code)
    spGenDiv <- rep(0, length(spAbundHarmMean))
    
    return(new('localComm',
               indSpecies = indSpecies,
               indTrait = indTrait,
               indSeqs = indSeqs,
               spGenDiv = spGenDiv,
               spAbund = spAbund,
               spTrait = spTrait,
               spAbundHarmMean = spAbundHarmMean,
               spLastOriginStep = spLastOriginStep,
               spExtinctionStep = spExtinctionStep))
}

#' @title The metacommunity of a `roleData` 
#' 
#' @description An S4 class to specify the state of a metacommunity
#'
#' @slot spAbund a numeric vector of the relative abundances for each species in 
#'     the meta community
#' @slot spTrait a numeric vector of the average trait value for each species in 
#'     the meta community

#' @rdname metaComm
#' @import methods
#' @export

setClass('metaComm',
         slots = c(spAbund = 'numeric',
                   spTrait = 'numeric'))


# constructor 
#' @rdname metaComm
#' @param spAbund spAbund
#' @param spTrait spTrait
#' @import methods
#' @export

metaComm <- function(spAbund, spTrait) {
    return(new('metaComm',
               spAbund = spAbund,
               spTrait = spTrait))
}
