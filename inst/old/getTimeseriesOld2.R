
#' @title Get a timeseries from a roleModel or roleExperiment as a data.frame 
#'
#' @param x a roleModel or roleExperiment object to get a stat from 
#' @param runNum an integer specifying the run number if a roleExperiment object is used
#' @param type a string specifying the type of values to get, either "summary_stats", or "model_values" 

#' @return a dataframe containing the timeseries 
#'
#' @export

setGeneric('getTimeseries', function(x, ...) standardGeneric('getTimeseries'), signature = 'x')
setMethod("getTimeseries", signature(x="roleModel"),
          function(x,type) {
            
            # get number of timeseries steps
            nseries <- length(x@timeseries)
            
            #create data frame with 0 rows and 5 columns
            ts_df <- data.frame(matrix(ncol = 0, nrow = nseries))
            
            #set row names to be the number of iterations as in 0, niterTimestep * 1, niterTimestep * 2 ... niterTimestep * nseries
            rownames(ts_df) <- x@arguments@niterTimestep * c(0:nseries)
            
            if(type == "summary_stats")
            {
              # get unique entropies and hill number types
              entropies <- unique(x@timeseries[[1]]@stats$entropy)
              hill_types <- unique(x@timeseries[[1]]@stats$type)
              
              # for every entropy and type
              all_stats_ts <- list()
              for(e in unique_e){
                for(t in unique_t){
                  
                  # for every step in the series
                  stat_ts <- c()
                  for(i in 1:nseries)
                  {
                    # get stats at the series 
                    stats <- x@timeseries[[i]]@stats
                    # add that stat series step
                    stat_ts <- c(stat_ts,subset(stats, type == t & entropy == e))
                  }
                  
                  # add stat ts to whole ts df
                  ts_df[,paste0(e,"_",t)] = stat_ts
                }
              }
            }
            
            else if(type == "model_values")
            {
              # due to how it's organized, the best way to access by name is by coercing 
              # everything into a named list then getting the name
              
              all_values_ts <- list()
              
              for(i in 1:nseries)
              {
                object_slots <- c(rep("localComm",6),rep("metaComm",2),rep("phylo",4))
                value_slots <- c("abundanceIndv","speciesIDsIndv","traitsIndv","abundanceSp","traitsSp","gDiversitiesSp",
                                 "abundanceSpM","metaTraitsSpM",
                                 "ntips","edges","lengths","alive")
                for(s in 1:length(value_slots)){
                  value_slot_name = value_slots[s]
                  all_values_ts[[value_slot_name]][[i]] <- slot(x@timeseries[[i]]@localComm, value_slot_name)
                }
              }
              
              for(s in 1:length(value_slots)){
                value_slot_name = value_slots[s]
                ts_df[,value_slot_name] <- I(all_values_ts[[value]])
              }
            }
            
            return(ts_df)
          }
)

setMethod("getTimeseries", signature(x="roleExperiment"),
          function(x,type,runNum) {
            
            model <- x@modelRuns[[runNum]]
            
            return(getTimeseries(model,type))
          }
)
