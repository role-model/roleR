#' @title The local community of a `roleData` 
#' 
#' @description An S4 class to specify the state of a local community
#' @slot indSpecies a numeric vector of the species IDs for each individual
#' For example indSpecies[2] is the ID of the species of individual 2
#' @slot indTrait a numeric vector of the trait values for each individual
#' @slot indSeqs a character vector of the gene sequences for each individual
#' @slot spGenDiv a numeric vector of the genetic diversities for each species
#' @slot spTrait a numeric vector of the mean trait value of each species
#' @slot spAbund a numeric vector of the abundance (number of individuals) of each species
#' @slot spAbundHarmMean numeric vector of the harmonic mean of species abundances
#' Equivalent to n / the reciprocal sum of abundances since last emergence where n is the number of steps in the reciprocal sum
#' Used in genetic simulation in roleR
#' @slot spLastOriginStep numeric vector of the last origin step of each species
#' The last origin step is the last step where a new species appeared in the local community (either totally new or after an extinction in the local)
#' Used in genetic simulation in roleR
#' @slot spExtinctionStep numeric vector of the last extinction (extirpation) step of each species
#' Used in genetic simulation in roleR
#' @slot equilibProp numeric something
#' @rdname localComm
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
#' @export

localComm <- function(indSpecies, indTrait, indSeqs, spGenDiv) { #indSppTrt, indSeqs, sppGenDiv
    
    # create the species indexed vectors from individual index vectors
    spAbund <- as.vector(tabulate(indSpecies))
    spTrait <- as.vector(tapply(indTrait, indSpecies, mean))
    
    # initialize other vectors that will be filled as the model runs
    spAbundHarmMean <- as.numeric(c())
    spLastOriginStep <- as.numeric(c())
    spExtinctionStep <- as.numeric(c())
    
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
#' @slot spAbund a numeric vector of the relative abundances for each species in the meta community
#' For example spAbundTrt[3] is the relative abundance of species 3
#' @slot spTrait a numeric vector of the average trait value for each species in the meta community

#' @rdname metaComm
#' @export

setClass('metaComm',
         slots = c(spAbund='numeric',
                   spTrait='numeric'))


# constructor 
#' @rdname metaComm
#' @param spAbund spAbund
#' @param spTrait spTrait
#' @export

metaComm <- function(spAbund,spTrait) {
    return(new('metaComm',
               spAbund=spAbund,
               spTrait=spTrait))
}