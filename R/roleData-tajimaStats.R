#' @title Summary statistics based on Tajima's D for `roleData` objects.
#' @description These functions summarize Tajima's D of the species in the 
#'     local community, and how it relates to nucleotide diversity (pi) and 
#'     abundance. `meanTajD` calculates the mean of Tajima's D across species; 
#'     `hillTajD` calculates trait-based Hill numbers treating Tajima's D as 
#'     a trait; and the `cor*` functions calculate Spearman correlations across
#'     species between Tajima's D and pi (`corTajDPi`), Tajima's D and 
#'     abundance (`corTajDAbund`), and abundance and pi (`corAbundPi`).
#' @param x The `roleData` object to calculate the statistics from
#' @param q The Hill number exponents. Can be a single value or a vector of integers.  
#' @param ... additional parameters, ignored
#' 
#' @details Tajima's D and pi come from the genetic simulation, so they are 
#'     `NA` for all species in a model that has not been run, and for species
#'     that are not in the local community at a snapshot. Species with `NA` 
#'     values are left out of every statistic. Correlations are `NA` when fewer 
#'     than 3 species remain or when either variable is constant, and `hillTajD`
#'     is `NA` when fewer than 2 species remain.
#'     
#'     For `hillTajD` species abundance is used as the weight and Tajima's D as 
#'     the trait, see `hillTrait` for the calculation.
#' 
#' @return A single number, or for `hillTajD` a vector with one number for each
#'     value of `q`
#' 
#' @include roleData.R
#' @importFrom stats cor
#' @rdname tajD-sumStats
#' @export

# mean tajima D
setGeneric('meanTajD', 
           def = function(x, ...) standardGeneric('meanTajD'), 
           signature = 'x')



#' meanTajD on roleData
#' @name meanTajD
#' @aliases meanTajD,roleData-method
#' @docType methods
#' @rdname tajD-sumStats

setMethod('meanTajD', 
          signature = 'roleData', 
          definition = function(x) {
              D <- x@localComm@spTajD
              D <- D[is.finite(D)]
              
              return(mean(D, na.rm = TRUE))
          }
)


# tajima D hill stats
#' @rdname tajD-sumStats
#' @export

setGeneric('hillTajD', 
           def = function(x, q = 1:4, ...) standardGeneric('hillTajD'), 
           signature = 'x')



#' hillTajD on roleData
#' @name hillTajD
#' @aliases hillTajD,roleData-method
#' @docType methods
#' @rdname tajD-sumStats

setMethod('hillTajD', 
          signature = 'roleData', 
          definition = function(x, q = 1:4) {
              X <- cbind(x@localComm@spAbund, x@localComm@spTajD)
              X <- X[!is.na(X[, 2]) & is.finite(X[, 2]) & X[, 1] > 0, , 
                     drop = FALSE]
              
              if(nrow(X) < 2) return(rep(NA_real_, length(q)))
              
              return(.hillDivTrait(X, q))
          }
)


# correlation between tajima D and pi
#' @rdname tajD-sumStats
#' @export

setGeneric('corTajDPi', 
           def = function(x, ...) standardGeneric('corTajDPi'), 
           signature = 'x')



#' corTajDPi on roleData
#' @name corTajDPi
#' @aliases corTajDPi,roleData-method
#' @docType methods
#' @rdname tajD-sumStats

setMethod('corTajDPi', 
          signature = 'roleData', 
          definition = function(x) {
              return(.rankCorrSpp(x@localComm@spTajD, x@localComm@spGenDiv))
          }
)


# correlation between tajima D and abundance
#' @rdname tajD-sumStats
#' @export

setGeneric('corTajDAbund', 
           def = function(x, ...) standardGeneric('corTajDAbund'), 
           signature = 'x')



#' corTajDAbund on roleData
#' @name corTajDAbund
#' @aliases corTajDAbund,roleData-method
#' @docType methods
#' @rdname tajD-sumStats

setMethod('corTajDAbund', 
          signature = 'roleData', 
          definition = function(x) {
              return(.rankCorrSpp(x@localComm@spTajD, x@localComm@spAbund))
          }
)


# correlation between abundance and pi
#' @rdname tajD-sumStats
#' @export

setGeneric('corAbundPi', 
           def = function(x, ...) standardGeneric('corAbundPi'), 
           signature = 'x')



#' corAbundPi on roleData
#' @name corAbundPi
#' @aliases corAbundPi,roleData-method
#' @docType methods
#' @rdname tajD-sumStats

setMethod('corAbundPi', 
          signature = 'roleData', 
          definition = function(x) {
              return(.rankCorrSpp(x@localComm@spAbund, x@localComm@spGenDiv))
          }
)


# Spearman correlation of two species-level vectors
# species with `NA` in either vector are dropped
# @param a,b numeric vectors of equal length with one element per species

.rankCorrSpp <- function(a, b) {
    ok <- !is.na(a) & !is.na(b) & is.finite(a) & is.finite(b)
    
    if(sum(ok) < 3) return(NA_real_)
    
    # a constant vector has no rank correlation
    if(length(unique(a[ok])) < 2 || length(unique(b[ok])) < 2) {
        return(NA_real_)
    }
    
    return(cor(a[ok], b[ok], method = "spearman"))
}
