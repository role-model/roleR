
#' @title Get a timeseries of a summary statistic as a vector or matrix
#'
#' @param x a roleModel or roleExperiment object to get a stat from 
#' @param runNum an integer specifying the run number if a roleExperiment object is used
#' @param type a string specifying the type of value to get, either "summary_stat", or "model_value" 
#' @param valueName a string specfiying the field to get i.e. abundanceIndv 
#' @param entropyN a integer specifying the Hill entropy to get a value for

#' @return a vector or matrix containing the timeseries 
#'
#' @export

x <- readRDS("test.roleexperiment")
x <- x@modelRuns[[1]]
type <- "model_value"
valueName <- "abundanceIndv"

# sets a parameter of roleParams
setGeneric('getTimeseries', function(x, ...) standardGeneric('getTimeseries'), signature = 'x')
setMethod("getTimeseries", signature(x="roleModel"),
          function(x,type,valueName,entropyN) {
            
            # get number of iterations of model
            niter <- length(x@paramValues$individuals_local)
            if(type == "summary_stat")
            {
              all_values_ts <- list()
              for(i in 1:niter)
              {
                #i = 1
                stat <- subset(x@timeseries@stats, hillType == name & hillEntropy == entropyN)
                stat <- stat[1,1]
                all_values_ts <- c(all_values_ts,stat)
              }
              return(all_values_ts)
            }
            
            else if(type == "model_value")
            {
              all_values_ts <- list()
              
              # due to how it's organized, the best way to access by name is by coercing 
              # everything into a named list then getting the name
              
              # 1:n_timesteps
              for(i in 1:10)
              {
                print(i)
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

              list_lists <- all_values_ts[["abundanceSp"]]
              #mat <- do.call("cbind",list_lists)
              #mat
              return(all_values_ts[[valueName]])
            }
          }
)

setGeneric('getTimeseries', function(x, ...) standardGeneric('getTimeseries'), signature = 'x')
setMethod("getTimeseries", signature(x="roleExperiment"),
          function(x,runNum,type,valueName,entropyN) {
            
            model <- x@modelRuns[[runNum]]
            
            return(getTimeseries(model,type,valueName,entropyN))
)
