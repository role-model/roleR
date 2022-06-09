
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

#x <- readRDS("test.roleexperiment")
#x <- x@modelRuns[[1]]
#type <- "model_value"
#valueName <- "abundanceIndv"

# sets a parameter of roleParams
setGeneric('getTimeseries', function(x, ...) standardGeneric('getTimeseries'), signature = 'x')
setMethod("getTimeseries", signature(x="roleModel"),
          function(x,type,valueName=NULL) {
            
            # get number of steps in timeseries
            nseries <- length(x@timeseries)
            all_values_ts <- list()
            
            # 
            if(type == "summary_stat")
            {
              unique_e <- unique(x@timeseries[[1]]@stats$entropy)
              unique_t <- unique(x@timeseries[[1]]@stats$type)
              all_stats_ts <- list()
              for(e in unique_e){
                for(t in unique_t){
                  stat_ts <- c()
                  for(i in 1:nseries)
                  {
                    stats <- x@timeseries[[i]]@stats
                    all_values_ts[[paste0(t,"_",e)]][[i]] <- subset(stats, type == t & entropy == e)[1,3]
                  }
                }
              }
              if(is.null(valueName)){
                all_values_ts <- do.call(rbind.data.frame, all_values_ts)
                all_values_ts <- t(all_values_ts)
                rownames(all_values_ts) <- NULL
                all_values_ts <- as.data.frame(all_values_ts)
                all_values_ts <- data.frame(sapply(all_values_ts, function(x) as.numeric(as.character(x))))
                return(all_values_ts)
              }
              else{
                return(as.numeric(all_values_ts[[valueName]]))
              }
            }
            
            else if(type == "model_value")
            {
              all_values_ts <- list()
              
              # due to how it's organized, the best way to access by name is by coercing 
              # everything into a named list then getting the name
              
              # 1:n_timesteps
              for(i in 1:nseries)
              {
                #print(i)
                all_values_ts[["abundanceIndv"]][[i]] <- x@timeseries[[i]]@localComm@abundanceIndv 
                all_values_ts[["speciesIDsIndv"]][[i]] <- x@timeseries[[i]]@localComm@speciesIDsIndv
                all_values_ts[["traitsIndv"]][[i]] <- x@timeseries[[i]]@localComm@traitsIndv
                all_values_ts[["abundanceSp"]][[i]]  <- x@timeseries[[i]]@localComm@abundanceSp
                all_values_ts[["traitsSp"]][[i]] <- x@timeseries[[i]]@localComm@traitsSp
                #all_values_ts[["sequencesSp"]][[i]] <- x@timeseries[[i]]@localComm@sequencesSp
                
                all_values_ts[["metaAbundanceSp"]][[i]] <- x@timeseries[[i]]@metaComm@abundanceSp
                all_values_ts[["metaTraitsSp"]][[i]] <- x@timeseries[[i]]@metaComm@traitsSp
                
                all_values_ts[["ntips"]][[i]] <- x@timeseries[[i]]@phylo@ntips
                all_values_ts[["edges"]][[i]] <- x@timeseries[[i]]@phylo@ntips
                all_values_ts[["lengths"]][[i]] <- x@timeseries[[i]]@phylo@ntips
                all_values_ts[["alive"]][[i]] <- x@timeseries[[i]]@phylo@ntips
              }
              
              if(is.null(valueName)){
                return(all_values_ts)
              }
              else{
                return(all_values_ts[[valueName]])
              }
            }
            
          }
)

setMethod("getTimeseries", signature(x="roleExperiment"),
          function(x,runNum,type,valueName=NULL) {
            
            model <- x@modelRuns[[runNum]]
            
            return(getTimeseries(model,type,valueName,entropyN))
          }
)
