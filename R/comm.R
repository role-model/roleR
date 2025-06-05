#' @title The local community of a `roleData` object
#' 
#' @description An S4 class to specify the state of a local community
#' @slot indSpecies a numeric vector of the species IDs for each individual
#' @slot indTrait a numeric vector of the trait values for each individual
#' @slot indSeqs a character vector of the gene sequences for each individual
#' @slot spAbund a numeric vector of abundances for each species 
#' @slot spGenDiv a numeric vector of the genetic diversities for each species
#' @slot spTajD a numeric vector of the Tajima D stat for each species
#' @slot spAbundHarmMean numeric vector of the harmonic mean of species 
#'    abundances
#' @slot spAbundHarmMean numeric vector of the harmonic mean of species 
#'    abundances
#' @slot spLastOriginStep numeric vector holding time of most recent origin of 
#'    each species in the local community
#' @slot gens numeric scalar of (fractional) number of generations transpired
#' @slot equilibProp numeric proportion of equilibrium achieved 
#' 
#' @rdname localComm
#' 
#' @import methods
#' 
#' @export

setClass('localComm',
         slots = c(indSpecies = 'numeric',
                   indTrait = 'matrix', 
                   indSeqs = 'character',
                   spAbund = 'numeric',
                   spGenDiv = 'numeric',
                   spTajD = 'numeric',
                   spAbundHarmMean = 'numeric',
                   spPropHarmMean = 'numeric',
                   spLastOriginStep = 'numeric',
                   gens = 'numeric',
                   equilibProp = 'numeric'))


# constructor 
#' @rdname localComm
#' @param indSpecies indSpecies
#' @param indTrait indTrait
#' @param indSeqs indSeqs
#' @import methods
#' @export

localComm <- function(indSpecies, indTrait, indSeqs) { 
    
    # create the species indexed vectors from individual index vectors 
    # and initialize with NA for missing data where appropriate
    
    spAbundHarmMean <- spAbund <- as.vector(tabulate(indSpecies))
    spPropHarmMean <- spAbundHarmMean / length(indSpecies)
    
    spLastOriginStep <- numeric(length(spAbundHarmMean))
    spLastOriginStep[spAbundHarmMean == 0] <- NA # if abund = 0 not there yet
    
    spGenDiv <- spTajD <- as.numeric(rep(NA, length(spAbundHarmMean)))
    
    return(new('localComm',
               indSpecies = indSpecies,
               indTrait = indTrait,
               indSeqs = indSeqs,
               spAbund = spAbund,
               spGenDiv = spGenDiv,
               spTajD = spTajD,
               spAbundHarmMean = spAbundHarmMean,
               spPropHarmMean = spPropHarmMean,
               spLastOriginStep = spLastOriginStep, 
               gens = 0, 
               equilibProp = 0))
}

#' @title The metacommunity of a `roleData` object
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
                   spTrait = 'matrix'))


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
