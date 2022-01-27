#' @title An S4 class to specify parameters of the RoLE model
#'
#' @slot params a named list of numerical values of the parameters
#' @slot timeseries a named list containing vectors of parameters as they vary across the simulation
#' Instead of specifying a prior distribution for a parameter, sample a parameter many times from a prior distribution of your choosing
#' and place it in the timeseries slot, where values are drawn from sequentially during the simulation 
#' i.e. sampling of parameters for the whole model run is done beforehand, due to how much easier and more flexible sampling
#' is in R as opposed to C++ 

#' @export

roleParams <- setClass('roleParams',

         slots = c(
           nsim = "numeric", 
           niter = "numeric",
           niter_timestep = "numeric",
           params = "list"
         ), 
         
         validity=function(object)
         {
           return(TRUE)
         })

# method defined as a constructor for roleParams
roleParams <- function(nsim, niter, niter_timestep = NULL) {
  params <- vector(mode = "list", length = nsim)
  for(i in 1:nsim)
  {
    p <- data.frame(matrix(ncol = 10, nrow = niter))
    colnames(p) <- c("speciation_meta", "extinction_meta", "trait_sigma", "species_meta", "individuals_meta",
                      "individuals_local", "dispersal_prob", "speciation_local")
    params[[i]] <- p
  }
  return(new("roleParams", nsim=nsim, niter=niter,params=params,niter_timestep=niter_timestep))
}

# sets a parameter of roleParams
setMethod("setParam", signature(x="roleParams"),
          function(x,param_name,value,sim_num) {
            
            # if no sim_num provided, apply to all sims
            # bad code but don't understand how to better do this in S4 
            if(is.null(sim_num))
            {
              for(i in 1:x@nsim)
              {
                if(length(value) == 1){ # if only one value specified, set all iterations
                  x@params[i]$param_name <- rep(value,x@niter)
                }
                else {
                  x@params[i]$param_name <- value
                }  
              }
            }
            
            # else apply only to the target sim_num
            else
            {
              if(length(value) == 1){ # if only one value specified, set all iterations
                x@params[[sim_num]][,param_name] <- rep(value,x@niter)
              }
              else {
                x@params[[sim_num]][,param_name] <- value
              }  
            }
          }
)

# method adding default params
setGeneric('addDefaultParams', function(x, ...) standardGeneric('addDefaultParams'), signature = 'x')
setMethod("addDefaultParams", signature(x="roleParams"),
          function(x) {
            #setParam(x, )
            
          }
)

# sets a parameter of roleParams
setGeneric('toCpp', function(x) standardGeneric('toCpp'), signature = 'x')
setMethod("toCpp", signature(x="roleParams"),
          function(x) {
            
          }
)