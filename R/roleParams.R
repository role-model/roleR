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
           nruns = "numeric", 
           niter = "numeric",
           values = "array",
           niter_timestep = "numeric"
         ), 
         
         validity=function(object)
         {
           return(TRUE)
         })

# method defined as a constructor for roleParams
roleParams <- function(nruns, niter, niter_timestep = NULL) {

  values <- array(list(),nruns)
  
  # set default zeroes
  for(i in 1:nruns)
  {
    values[[i]][["individuals_local"]] <- c(0)
    values[[i]][["individuals_meta"]] <- c(0)
    values[[i]][["species_meta"]] <- c(0)
    values[[i]][["speciation_local"]] <- c(0)
    values[[i]][["speciation_meta"]] <- c(0)
    values[[i]][["trait_sigma"]] <- c(0)
    values[[i]][["env_sigma"]] <- c(0)
    values[[i]][["comp_sigma"]] <- c(0)
    values[[i]][["dispersal_prob"]] <- c(0)
    values[[i]][["mutation_rate"]] <- c(0)
    values[[i]][["equilib_escape"]] <- c(0)
  }

    #names(values[[1]]) <- c("individuals_local", "individuals_meta", "species_meta", 
    #              "speciation_local", "speciation_meta", 
     #             "trait_sigma", "env_sigma", "comp_sigma", 
     #             "dispersal_prob", "mutation_rate", "equilib_escape")
  return(new("roleParams", nsim=nsim, niter=niter,values=values,niter_timestep=niter_timestep))
}

# sets a parameter of roleParams
setMethod("setParam", signature(x="roleParams"),
          function(x,paramName,values,runNum) {
            
            # if no runNum provided, apply to all sims
            # bad code but don't understand how to better do this in S4 
            if(is.null(sim_num))
            {
              for(i in 1:x@nsim)
              {
                if(length(value) == 1){ # if only one value specified, set all iterations to that one value
                  x@values[[i]][[param_name]] <- rep(values,x@niter)
                }
                else { # otherwise set all iterations to all values
                  x@values[[i]][[param_name]] <- values
                }  
              }
            }
            
            # else apply only to the target sim_num
            else
            {
              if(length(value) == 1){ # if only one value specified, set all iterations to that one value
                x@values[[i]][[param_name]] <- rep(values,x@niter)
              }
              else { # otherwise set all iterations to all values
                x@values[[i]][[param_name]] <- values
              } 
            }
          }
)

# method adding default params
setGeneric('setDefaultParams', function(x, ...) standardGeneric('setDefaultParams'), signature = 'x')
setMethod("setDefaultParams", signature(x="roleParams"),
          function(x) {
            setParam(x,"species_meta",15)
            setParam(x,"individuals_meta",500)
            setParam(x,"individuals_local",100)
            setParam(x,"dispersal_prob",0.5)
            setParam(x,"speciation_local",0.1)
            setParam(x,"extinction_meta",0.1)
            setParam(x,"speciation_meta",0.1)
            setParam(x,"trait_sigma",100)
            setParam(x,"sigma_e",0.1)
            setParam(x,"sigma_c",0.1)
            setParam(x,"trait_z",0.5)
          }
)

# sets a parameter of roleParams
setGeneric('toCpp', function(x) standardGeneric('toCpp'), signature = 'x')
setMethod("toCpp", signature(x="roleParams"),
          function(x) {
            
          }
)