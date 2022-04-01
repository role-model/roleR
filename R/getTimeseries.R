
#' @title Get a timeseries of a summary statistic as a vector or matrix
#'
#' @param x a roleModel or roleExperiment object to get a stat from 
#' @param runNum an integer specifying the run number if a roleExperiment object is used
#' @param type a string specifying the type of value to get, either "summary_stat", or "model_value" 
#' @param entropyN a integer specifying the Hill entropy to get a value for

#' @return a vector or matrix containing the timeseries 
#'
#' @export

# sets a parameter of roleParams
setGeneric('getTimeseries', function(x, ...) standardGeneric('getTimeseries'), signature = 'x')
setMethod("getTimeseries", signature(x="roleModel"),
          function(x,type,valueName,entropyN) {
            
            x <- dummyModel(R=FALSE,run=TRUE)
            
            type <- "model_value"
            valueName <- "abundance_local"
              
            if(type == "summary_stat")
            {
              niter <- length(x@paramValues[["speciation_local"]])
              out <- list()
              for(i in 1:niter)
              {
                stat <- subset(x@timeseries@stats, hillType == name & hillEntropy == entropyN)
                stat <- stat[1,1]
                out <- c(out,stat)
              }
              return(out)
            }
            
            else if(type == "model_value")
            {
              #niter <- x@paramValues[["niter"]]
            
              niter <- length(x@paramValues[["speciation_local"]])
              niter <- 1
              x@paramValues[["speciation_local"]]
              all_values_ts <- list()
              
              # due to how it's organized, the best way to access by name is by coercing 
              # everything into a named list then getting the name
              View(x@timeseries[1][[1]])
              for(i in 1:niter)
              {
                all_values_ts[["abundanceIndv"]][[i]] <- x@timeseries[i][[1]]@localComm@abundanceIndv 
                all_values_ts[["speciesIDsIndv"]][[i]] <- x@timeseries[i][[1]]@localComm@speciesIDsIndv
                all_values_ts[["traitsIndv"]][[i]] <- x@timeseries[i][[1]]@localComm@traitsIndv
                all_values_ts[["abundanceSp"]][[i]]  <- x@timeseries[i][[1]]@localComm@abundanceSp
                all_values_ts[["traitsSp"]][[i]] <- x@timeseries[i][[1]]@localComm@traitsSp
                #all_values_ts[["sequencesSp"]][[i]] <- x@timeseries[i][[1]]@localComm@sequencesSp
                
                all_values_ts[["metaAbundanceSp"]][[i]] <- x@timeseries[i][[1]]@metaComm@abundanceSp
                all_values_ts[["metaTraitsSp"]][[i]] <- x@timeseries[i][[1]]@metaComm@traitsSp
                
                all_values_ts[["ntips"]][[i]] <- x@timeseries[i][[1]]@phylo@ntips
                all_values_ts[["edges"]][[i]] <- x@timeseries[i][[1]]@phylo@ntips
                all_values_ts[["lengths"]][[i]] <- x@timeseries[i][[1]]@phylo@ntips
                all_values_ts[["alive"]][[i]] <- x@timeseries[i][[1]]@phylo@ntips
                
                #all_values <- rbind(x@timeseries[i])
              }
              all_values_ts
              all_values_ts[["abundanceIndv"]]
              
              return(all_values_ts[[valueName]])
            }
          }
)
