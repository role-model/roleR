#' @title An S4 class to specify parameters of the RoLE model
#'
#' @slot params a named list of numerical values of the parameters
#' @slot timeseries a named list containing vectors of parameters as they vary across the simulation
#' Instead of specifying a prior distribution for a parameter, sample a parameter many times from a prior distribution of your choosing
#' and place it in the timeseries slot, where values are drawn from sequentially during the simulation 
#' i.e. sampling of parameters for the whole model run is done beforehand, due to how much easier and more flexible sampling
#' is in R as opposed to C++ 

#' @export

roleParameters <- setClass('roleParameters',
         slots = c(params = 'list',
                   timeseries = 'list'),

         prototype=list(
           params = list(species_meta = 15,
                            individuals_meta = 500,
                            individuals_local = 100,
                            dispersal_prob = 0.5,
                            speciation_local = 0.1,
                            extinction_meta = 0.8,
                            speciation_meta = 1,
                            trait_sigma = 0.1),
           timeseries = list(species_meta = NULL,
                            individuals_meta = NULL,
                            individuals_local = NULL,
                            dispersal_prob = NULL,
                            speciation_local = NULL,
                            extinction_meta = NULL,
                            speciation_meta = NULL,
                            trait_sigma = NULL)
          ),
         
         validity=function(object)
         {
           return(TRUE)
         })