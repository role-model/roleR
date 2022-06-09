
#' @title Get a timeseries from a roleModel or roleExperiment as a dataframe 
#'
#' @param x a roleModel or roleExperiment object to get a timeseries from
#' @param fun a function that applies on roleData objects to gather data 
#' apply applies the fun across the list of roleData in a roleModel
#' @param lastColAsTitle specify to use the last column gathered by the fun as column names for the output dataframe

#' @return a dataframe containing the timeseries 
#'
#' @export

setGeneric('getTimeseries', function(x, ...) standardGeneric('getTimeseries'), signature = 'x')
setMethod("getTimeseries", signature(x="roleModel"),
          function(x,fun,lastColAsTitle) {
            
            ts_lists <- lapply(x@timeseries,fun)
            ts_df <-  t(as.data.frame(do.call(cbind, ts_lists)))
            
            if(lastColAsTitle){
              names <- unlist(ts_df[1,ncol(ts_df)])
              ts_df <- ts_df[,1:(ncol(ts_df) - 1)]
              colnames(ts_df) <- unlist(ts_df[1,ncol(ts_df)])
            }
            
            return(ts_df) 
          }
)

role <- readRDS("data/example_out_in/test.roleexperiment")
model <- role@modelRuns[[1]]

# get one value
getAbundanceSp <- function(roledata){
  return(list(roledata@iterNum,roledata@localComm@abundanceSp))
} 

# get all abundances and name them with last col
getAbundances <- function(roledata){
  return(list(roledata@iterNum,roledata@localComm@abundanceSp,roledata@metaComm@abundanceSp,c("iterNum","abundanceLocal","abundanceMeta")))
} 

# get all abundances and name them with last col
getTraits <- function(roledata){
  return(list(roledata@iterNum,roledata@localComm@traitsIndv,roledata@localComm@traitsSp,roledata@metaComm@traitsSpM,
              c("iterNum","traitsIndv","traitsSp","traitsSpM")))
} 

# get all abundances and name them with last col
getTraitMean <- function(roledata){
  return(list(roledata@iterNum,mean(roledata@localComm@traitsIndv),
              c("iterNum","traitsIndvMean")))
} 

setMethod("getTimeseries", signature(x="roleExperiment"),
          function(x,fun,runNum=NULL) {
            
            # if supplied runNum
            if(!is.null(runNum)){
              model <- x@modelRuns[[runNum]]
              return(getTimeseries(model,fun))
            }
            
            # else return list of all model timeseries
            all_models_ts <- list()
            for(model in x@modelRuns)
            {
              all_models_ts <- list(all_models_ts,getTimeseries(model,fun))
            }
          }
)
