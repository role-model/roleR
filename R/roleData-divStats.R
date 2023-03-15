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

#' @rdname div-sumStats
#' @export

# abundance hill stats
#' @rdname div-sumStats
#' @export

setGeneric('hillAbund', 
           def = function(x, q = 1:4, ...) standardGeneric('hillAbund'), 
           signature = 'x')


setMethod('hillAbund', 
          signature = 'roleData', 
          definition = function(x, q = 1:4) {
              X <- rawAbundance(x)[[1]]
              X <- X[X > 0]
              
              return(.hillDivSAD(X, q))
          }
)


# function for abundance-based
#' @param n is a vector of species abundances

.hillDivSAD <- function(n, q) {
    n <- n/sum(n)
    
    hill <- outer(n, q, '^')
    hill <- colSums(hill)^(1 / (1 - q))
    hill[q == 1] <- exp(sum(-n * log(n)))
    
    return(hill)
}



# total abundance
#' @rdname total-abundance
#' @export

setGeneric('totalN', 
           def = function(x, ...) standardGeneric('totalN'), 
           signature = 'x')


setMethod('totalN', 
          signature = 'roleData', 
          definition = function(x) {
              X <- rawAbundance(x)[[1]]
              return(sum(X))
          }
)
# function for total abundance
#' @param n is a vector of species abundances

.totalAbundSAD <- function(n) {
    
    return(sum(n))
    
}


# genetic hill stats
#' @rdname div-sumStats
#' @export

setGeneric('hillGenetic', 
           def = function(x, q = 1:4, ...) standardGeneric('hillGenetic'), 
           signature = 'x')


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


# function for trait-based
#' @param X a matrix with first column = abund, second column = traits
#' @param traits is a vector of traits

.hillDivTrait <- function(X, q) {
    n <- X[, 1]
    traits <- X[, 2]
    
    p <- n / sum(n)
    dij <- as.matrix(dist(traits))
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

setMethod('richness', 
          signature = 'roleData', 
          definition = function(x) {
              length(unique(x@localComm@indSpecies))
          }
)
