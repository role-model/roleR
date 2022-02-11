
#' @title Get the state of a model at a timestep as a roleData object
#'
#' @param x a roleModel or roleSim object to get data from
#' @param runNum an integer specifying the run number if a roleSim object is used
#' @param iterNum an integer specifying the model iteration to get data from - must be a multiple
#' of the niter_timestep parameter and cannot be greater than the niter parameter
#' @param iterAsPercentile an integer specifying the model iteration to get data from - must be a multiple
#' of the niter_timestep parameter and cannot be greater than the niter parameter

#' @examples
#' GetStateData(x=sim, runNum = 1, iter = 50)
#' GetStateData(x=model, runNum = 1, iter = 0.5, iterAsPercentile=TRUE)
#' 
#' @return a roleData object
#'
#' @export
#' NOTE - consider allowing users to specify a percentile to get a state from i.e. 0.5 gets state halfway through the iterations

setGeneric('getStateData', function(x, ...) standardGeneric('getStateData'), signature = 'x')
setMethod("getStateData", signature(x="roleSim"),
          function(x,runNum,iter,iterAsPercentile) {
            # iterNum / niter_timestep gets the index within the timeseries of the iterNum specified
            if(!iterAsPercentile){
              return(x@modelRuns[runNum]@timeseries[iter/x@params@niter_timestep])
            }
            else{
              # iter becomes a percentile of niter
              iter <- x@params@niter * iter
              return(x@modelRuns[runNum]@timeseries[iter/x@params@niter_timestep])
            }
          }
)
setMethod("getStateData", signature(x="roleModel"),
          function(x,iter) {
            # iterNum / niter_timestep gets the index within the timeseries of the iterNum specified
            return(x@timeseries[iter/x@params[["niter_timestep"]]])
          }
)