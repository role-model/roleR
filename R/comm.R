#' @title Local community
#' 
#' @description An S4 class to specify the state of a local community
#'
#' @slot indSppTrt numeric matrix with rows for individuals; column 1 is species
#'     ID (matches row number in `metaComm@sppData`); column 2 is trait value; 
#'     column 3 is sequence (maybe? maybe shouldn't have it)
#' @slot indSeqs character matrix with rows for individuals; column 1 is 
#'     sequences for each individual
#' @slot sppGenDiv a numeric matrix with rows for species; column 1 is species
#'     level genetic diversity
#' 
#' @rdname localComm
#' @export

setClass('localComm',
         slots = c(#indSppTrt = 'matrix',
                   indSpecies = 'numeric',
                   indTrait = 'numeric',
                   indSeqs = 'character',
                   spGenDiv = 'numeric',
                   #spAbundTrt = 'matrix',
                   spTrait = 'numeric',
                   spAbund = 'numeric'
                    ))


# constructor 
#' @rdname localComm
#' @export

localComm <- function(indSpecies, indTrait, indSeqs, spGenDiv) { #indSppTrt, indSeqs, sppGenDiv
    
    # create the species indexed vectors
    spAbund <- rep(NA, 10000)
    spAbund[1] <- sum(indSpecies == 1)
    spTrait <- rep(NA, 10000)
    spTrait[1] <- mean(indTrait)
    
    #spAbundTrt <- matrix(data=NA, 10000, 2);
    #spAbundTrt[1,1] <- length(indSppTrt)
    #spAbundTrt[1,2] <- mean(indSppTrt[,2])
    
    return(new('localComm',
             indSpecies = indSpecies,
             indTrait = indTrait,
             indSeqs = indSeqs,
             spGenDiv = spGenDiv,
             spAbund = spAbund,
             spTrait = spTrait))
}

#' @title Meta community
#' 
#' @description An S4 class to specify the state of a meta community
#'
#' @slot sppAbundTrt numeric matrix with rows for species (row number is 
#'     species ID); column 1 is relative abundance; column 2 is trait mean
#' 
#' @rdname metaComm
#' @export

setClass('metaComm',
         slots = c(#sppAbundTrt = 'matrix'),
                spAbund='numeric',
                spTrait='numeric'))


# constructor 
#' @rdname metaComm
#' @export

metaComm <- function(spAbund,spTrait) {
  return(new('metaComm',
             spAbund=spAbund,
             spTrait=spTrait))
}
#install.packages("Rcpp")
#library(Rcpp)
#Rcpp.package.skeleton(name = "roleR")

#Rcpp::Rcpp.package.skeleton(name="roleR2")
