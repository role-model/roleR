#' @title Add Hill statistics to all saved timesteps within a roleSim or roleModel
#'
#' @description Hill number diversities for species abundance, trait, population
#' genetic, and phylogenetic data
#'
#' @param x a roleSim/roleModel object 
#' @param entrScales a vector of Renyei entropy scales to use when calculating Hill numbers
#'
#' @examples
#' addHillStats(sim,c(2,3,4))
#' 
#' @export

setGeneric('addHillStats', function(x, ...) standardGeneric('addHillStats'), signature = 'x')
setMethod("addHillStats", signature(x="roleExperiment"),
          function(x,entrScales) {
            # for each run in sim runs
            for(r in 1:length(x@modelRuns))
            {
              #r = 1
              # for each data object in the timeseries
              for(d in 1:length(x@modelRuns[[r]]@timeseries))
              {
                #d = 1
                data <- x@modelRuns[[r]]@timeseries[[d]]
                
                # create df of 3 cols, first is type "trait" "phylo" "abundance"
                # second is entropy scale, third is value
                stats <- matrix(NA, ncol = 3)
                
                # abundance Hills
                hill <- computeHill(data@localComm@abundanceSp,type="abundance",entropies = entrScales)
                rows <- cbind(rep("abundance",3),entrScales,hill)
                stats <- rbind(stats, rows)
                
                # trait Hills
                hill <- computeHill(data@localComm@abundanceSp,data@localComm@traitsSp,type="trait",entropies = entrScales)
                rows <- cbind(rep("trait",3),entrScales,hill)
                stats <- rbind(stats, rows)
              
                # phylo Hills
                # hill <- computeHill(data@localComm@abundanceSp,data@localComm@phylo,type="phylo",entropies = entrScales)
                # rows <- cbind(rep("phylo",3),entropies,hill)
                # stats <- rbind(stats, rows)
                
                stats <- as.data.frame(stats)
                stats <- stats[-1,]
                colnames(stats) <- c("type","entropy","value")
                
                x@modelRuns[[r]]@timeseries[[d]]@stats <- stats
              }
            }
            
            return(x)
          }
)
setMethod("addHillStats", signature(x="roleModel"),
          function(x,entrScales) {
          }
)

# user-inaccessible helper function to compute Hill numbers
computeHill <- function(x,traits=NULL,phylo=NULL,type,entropies){
  if(type == "abundance")
  {
    # renormalize abundance from 0-1 
    x <- x/sum(x)
    # raise every element in abundance vector to every exponent (rows equal to len(x), cols equal to len(q))
    hill <- outer(x, entropies, '^')
    # sum cols and rise to exponent
    hill <- colSums(hill)^(1 / (1 - entropies))
    hill[hill == 1] <- exp(sum(-x * log(x)))
    
    return(hill)
  }
  else if(type == "trait")
  {
    p <- x / sum(x)
    # distances of traits
    dij <- as.matrix(dist(traits))
    # matrix multiply to get Q 
    Q <- as.vector(p %*% dij %*% p)
    # all the ps times all the ps - distance between them / Q 
    a <- outer(p, p, '*') / Q
    
    # exponentiate the Q matrix 
    # as Q approaches 1 it approaches the shannon entropy 
    Hk <- sapply(entropies, function(qk) {
      if(qk == 1) {
        return(exp(- sum(dij * a * log(a))))
      } else {
        return(sum(dij * a^qk)^(1 / (1 - qk)))
      }
    })
    
    # convert to generalized (effective number) Hill number 
    # i.e. normalizing by total number of species if every species was equally abundant
    # converts from observed to effective based on balance of dominance and rarity 
    D <- sqrt(Hk / Q)
    return(D)
  }
  else if(type == "phylo")
  {
  }
}