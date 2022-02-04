#' @title Get a timeseries of a summary statistic, population etc. from a
#' roleModel object
#'
#' @param model object of class \code{roleModel}
#' @param param
#' @param simNum
#'
#' @return an object
#'
#' @export

getTimeseries <- function(model, simrun, stat)
{
  # initialize out vector with length equal to the number of time steps
  out <- numeric(1:length(model@timeseries))
  # for each time series 
  for(i in 1:length(model$timeseries))
  {
    out[i] <- model$timeseries[i]$params$values[param]
  }
  return(out)
}

#' @title Write a ROLE sim object or ROLE model object to a serialized R file for transfer, sharing, 
#' or to use with Docker. File contents can only be read within R. 
#'
#' @param model object of class \code{roleModel}
#' @param dir the directory to write to i.e. "data/models" 
#' @param fileName the name of the new file
#' @param saveTxt specifies whether to save a .txt describing some model information
#' for viewing outside of R}
#'
#' @return an object of class \code{roleModel}
#'
#' @export

writeModel <- function(model, dir = NULL, fileName, saveTxt = TRUE)
{
  if(is.null(dir)){
    dir = getwd()
  }
  saveRDS(model,paste0(dir,"/",fileName,".role"))
  if(saveTxt){
    con <- file(paste0(dir, "/", fileName, "_info", ".txt"))
    writeLines(c("Info for ROLE model object","Created by Jacob Idec"), con)
    close(con)
  }
}

#' @title Hill number diversities
#'
#' @description Hill number diversities for species abundance, trait, population
#' genetic, and phylogenetic data
#'
#'
#' @param x a vector of species abundances
#' @param q a vector of Renyei entropy scales
#' @param traits a vector of species trait values, must be in the same order as
#' \code{x}, no internal checks are provided, so make sure this is correct before
#' running the function!
#'
#' @details Stub
#'
#' @examples
#' hillDivSAD(1:10, 0:3)
#' hillDivTrait(1:10, 1:10, 0:3)
#'
#' @return A vector, same length as \code{q}, giving diversities for each
#' value of \code{q}
#' 
#' @export

# method adding default params
setGeneric('addHillStats', function(x, ...) standardGeneric('addDefaultParams'), signature = 'x')
setMethod("addHillStats", signature(x="roleSim"),
          function(x,re_scales) {
            
          }
)
setMethod("addHillStats", signature(x="roleModel"),
          function(x,re_scales) {
            
            # for each run in sim runs
            for(r in 1:length(x@runs))
            {
              # create 
              stats <- numeric(re_scales*4)
              model <- runs[r]
              model
            }
          }
)

addHillStats <- function(sim, q) {
  x <- x/sum(x)
  
  hill <- outer(x, q, '^')
  hill <- colSums(hill)^(1 / (1 - q))
  hill[hill == 1] <- exp(sum(-x * log(x)))
  
  return(hill)
}


# test
# hillDivSAD(1:10, 0:3)
# sapply(0:3, function(qi) hillR::hill_taxa(1:10, qi))




#' @rdname hillDiv
#' @export
hillDivTrait <- function(x, traits, q) {
  p <- x / sum(x)
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