#' @title Local community
#' 
#' @description An S4 class to specify the state of a local community

#' @slot indSpecies numeric vector of the species IDs for each individual
#' For example indSpecies[1] is the ID of the species of individual 1
#' @slot indTrait numeric vector of the trait values for each individual
#' @slot indSeqs character vector of the gene sequences for each individual
#'     
#' @slot spGenDiv numeric vector of the genetic diversities for each species
#' @slot spTrait numeric vector of the mean trait value of each species
#' @slot spAbund numeric vector of the abundance (number of individuals) of each species
#' 
#' @slot spAbundHarmMean numeric vector of the harmonic mean of species abundances
#' Equivalent to n / the reciprocal sum of abundances since last emergence where n is the number of steps in the reciprocal sum
#' Used in genetic simulation in roleR
#' @slot spLastOriginStep numeric vector of the last origin step of each species
#' The last origin step is the last step where a new species appeared in the local community (either totally new or after local extinction)
#' Used in genetic simulation in roleR
#' @slot spExtinctionStep numeric vector of the last extinction step of each species
#' Used in genetic simulation in roleR
#' 
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
                   spExtinctionStep = 'numeric'
                    ))


# constructor 
#' @rdname localComm
#' @export

localComm <- function(indSpecies, indTrait, indSeqs, spGenDiv) { #indSppTrt, indSeqs, sppGenDiv
    
    # create the species indexed vectors from individual index vectors
    spAbund <- rep(0, 10000)
    temp <- tabulate(indSpecies)
    spAbund[1:length(temp)] <- temp
    
    spTrait <- rep(0, 10000)
    temp <- tapply(indTrait, indSpecies, mean)
    spTrait[as.numeric(names(temp))] <- as.numeric(temp)
    
    #spAbundTrt <- matrix(data=NA, 10000, 2);
    #spAbundTrt[1,1] <- length(indSppTrt)
    #spAbundTrt[1,2] <- mean(indSppTrt[,2])
    
    # todo - verify that this initial state is correct 
    spAbundHarmMean <- rep(0, 10000)
    spLastOriginStep <- rep(0, 10000)
    spExtinctionStep <- rep(0, 10000)
    
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

#' @title Meta community
#' 
#' @description An S4 class to specify the state of a meta community
#'
#' @slot spAbund numeric vector of relative abundances for each species in the meta community
#' For example spAbundTrt[1] is the relative abundance of species 1
#' @slot spTrait numeric vector of trait values for each species in the meta community

#' @rdname metaComm
#' @export

setClass('metaComm',
         slots = c(spAbund='numeric',
                spTrait='numeric'))


# constructor 
#' @rdname metaComm
#' @export

metaComm <- function(spAbund,spTrait) {
  return(new('metaComm',
             spAbund=spAbund,
             spTrait=spTrait))
}
