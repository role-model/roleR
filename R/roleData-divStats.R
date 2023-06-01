#' @title Hill number summary statistics for `roleData` objects.
#' @description These functions calculate Hill numbers. The naming convention 
#'     indicates what type of Hill number each function calculates; e.g. 
#'     `hillAbund` calculates abundance-based Hill numbers, while `hillTrait` 
#'     calculates trait-based Hill numbers.
#' @param x The `roleData` object to calculate Hill numbers or richness from
#' @param q The Hill number exponents. Can be a single value or a vector of integers.  
#' @param ... additional parameters, ignored
#' 
#' @details See Gaggiotti et al. 2018 for an overview of Hill numbers in the context of ecology and evolution. 
#' @references Gaggiotti, Oscar E, Anne Chao, Pedro Peres-Neto, Chun-Huo Chiu, Christine Edwards, Marie-Josée Fortin, Lou Jost, Christopher M Richards, and Kimberly A Selkoe. “Diversity from Genes to Ecosystems: A Unifying Framework to Study Variation across Biological Metrics and Scales.” Evol. Appl. 11, no. 7 (2018): 1176–93. https://doi.org/10.1111/eva.12593.
#' @include roleData.R
#' @rdname div-sumStats
#' @export


setGeneric('hillAbund', 
           def = function(x, q = 1:4, ...) standardGeneric('hillAbund'), 
           signature = 'x')


#' HillAbund on roleData
#' @name hillAbund
#' @aliases hillAbund,roleData-method
#' @docType methods
#' @rdname div-sumStats

setMethod('hillAbund', 
          signature = 'roleData', 
          definition = function(x, q = 1:4) {
              X <- rawAbundance(x)[[1]]
              X <- X[X > 0]
              
              return(.hillDivSAD(X, q))
          }
)


#' @title HillDivSAD
#' function for abundance-based
#'@param q  hill order
#'
#' @param n is a vector of species abundances

.hillDivSAD <- function(n, q) {
    # Remove zero values, primarily for genDiv where zero is a valid value
    n <- n[n != 0]
    n <- n/sum(n)
    
    hill <- outer(n, q, '^')
    hill <- colSums(hill)^(1 / (1 - q))
    hill[q == 1] <- exp(sum(-n * log(n)))
    
    return(hill)
}




# genetic hill stats
#' @rdname div-sumStats
#' @export

setGeneric('hillGenetic', 
           def = function(x, q = 1:4, ...) standardGeneric('hillGenetic'), 
           signature = 'x')



#' hillGenetic on roleData
#' @name hillGenetic
#' @aliases hillGenetic,roleData-method
#' @docType methods
#' @rdname div-sumStats

setMethod('hillGenetic', 
          signature = 'roleData', 
          definition = function(x, q = 1:4) {
              X <- rawGenDiv(x)[[1]]
              X <- X[!is.na(X)]
              
              return(.hillDivSAD(X, q))
          }
)


# trait hill stats
#' @rdname div-sumStats
#' @export

setGeneric('hillTrait', 
           def = function(x, q = 1:4, ...) standardGeneric('hillTrait'), 
           signature = 'x')



#' hillTrait on roleData
#' @name hillTrait
#' @aliases hillTrait,roleData-method
#' @docType methods
#' @rdname div-sumStats
setMethod('hillTrait', 
          signature = 'roleData', 
          definition = function(x, q = 1:4) {
              spp <- rawSppID(x)[[1]]
              trt <- rawTraits(x)[[1]]
              
              X <- lapply(unique(spp), function(i) {
                  xs <- sum(spp == i)
                  xt <- mean(trt[spp == i])
                  
                  return(c(xs, xt))
              })
              X <- do.call(rbind, X)
              
              return(.hillDivTrait(X, q))
          }
)

#' @title HillDivTrait
#' function for trait-based
#' @param X a matrix with first column = abund, second column = traits
#' @param q order for hill number

.hillDivTrait <- function(X, q) {
    
    n <- X[, 1]
    traits <- X[, 2]
    
    p <- n / sum(n)
    dij <- as.matrix(dist(traits))
    
    # catch case where species richness = 0
    if(nrow(X) == 1) {
        dij[1,1] <- 1
    }
    Q <- as.vector(p %*% dij %*% p)
    a <- outer(p, p, '*') / Q
    
    Hk <- sapply(q, function(qk) {
        if(qk == 1) {
            return(exp(- sum(dij * a * log(a))))
        } else {
            return(sum(dij * a^qk)^(1 / (1 - qk)))
        }
    })
    
    D <- sqrt(Hk / Q)
    
    
    return(D)
}


# phylo hill (place holder)
#' @rdname div-sumStats
#' @export

setGeneric('hillPhylo', 
           def = function(x, q = 1:4, ...) standardGeneric('hillPhylo'), 
           signature = 'x')


#' hillPhylo on roleData
#' @name hillPhylo
#' @aliases hillPhylo,roleData-method
#' @docType methods
#' @rdname div-sumStats
setMethod('hillPhylo', 
          signature = 'roleData', 
          definition = function(x, q = 1:4) {
              return(q)
          }
)



# species richness
#' @rdname div-sumStats
#' @export

setGeneric('richness', 
           def = function(x, ...) standardGeneric('richness'), 
           signature = 'x')


#' richness on roleData
#' @name richness
#' @aliases richness,roleData-method
#' @docType methods
#' @rdname div-sumStats
setMethod('richness', 
          signature = 'roleData', 
          definition = function(x) {
              length(unique(x@localComm@indSpecies))
          }
)
