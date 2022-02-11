
#' @title Get a timeseries of a summary statistic as a vector or matrix
#'
#' @param x a roleModel or roleSim object to get a stat from 
#' @param runNum an integer specifying the run number if a roleSim object is used
#' @param type a string specifying the type of value to get, either "summary_stat", or "model_value" 
#' @param entropyN a integer specifying the Hill entropy to get a value for

#' @return a vector or matrix containing the timeseries 
#'
#' @export

# sets a parameter of roleParams
setGeneric('getTimeseries', function(x, ...) standardGeneric('getTimeseries'), signature = 'x')
setMethod("getTimeseries", signature(x="roleModel"),
          function(x,type,valueName,entropyN) {
            
            if(type == "summary_stat")
            {
              niter <- x@params[["niter"]]
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
              niter <- x@params[["niter"]]
              
              all_values_ts <- list()
              
              # due to how it's organized, the best way to access by name is by coercing 
              # everything into a named list then getting the name
              for(i in 1:niter)
              {
                all_values_ts[["abundanceIndv"]][[i]] <- x@timeseries[i]@localComm@abundanceIndv 
                all_values_ts[["speciesIDsIndv"]][[i]] <- x@timeseries[i]@localComm@speciesIDsIndv
                all_values_ts[["traitsIndv"]][[i]] <- x@timeseries[i]@localComm@traitsIndv
                all_values_ts[["abundanceSp"]][[i]]  <- x@timeseries[i]@localComm@abundanceSp
                all_values_ts[["traitsSp"]][[i]] <- x@timeseries[i]@localComm@traitsSp
                all_values_ts[["sequencesSp"]][[i]] <- x@timeseries[i]@localComm@sequencesSp
                
                all_values_ts[["metaAbundanceSp"]][[i]] <- x@timeseries[i]@metaComm@abundanceSp
                all_values_ts[["metaTraitsSp"]][[i]] <- x@timeseries[i]@metaComm@traitsSp
                
                all_values_ts[["ntips"]][[i]] <- x@timeseries[i]@phylo@ntips
                all_values_ts[["edges"]][[i]] <- x@timeseries[i]@phylo@ntips
                all_values_ts[["lengths"]][[i]] <- x@timeseries[i]@phylo@ntips
                all_values_ts[["alive"]][[i]] <- x@timeseries[i]@phylo@ntips
                
                all_values <- rbind(x@timeseries[i])
              }
              
              return(all_values_ts[[valueName]])
            }
          }
)
