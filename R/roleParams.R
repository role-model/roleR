#' @title An S4 class to represent the parameters of one roleModel
#' @description has a slot for each parameter, each is a vector containing either one value or niter values 
#' 
#' @slot individuals_local
#' @slot individuals_meta
#' @slot species_meta
#' @slot speciation_local
#' @slot speciation_meta
#' @slot extinction_meta
#' @slot trait_sigma
#' @slot env_sigma
#' @slot comp_sigma
#' @slot dispersal_prob
#' @slot mutation_rate
#' @slot equilib_escape
#' @slot num_basepairs
#' @slot init_type a single character string either "oceanic_island" or "bridge_island"
#' @export

roleParams <- setClass('roleParams',
                          slots = c(
                            individuals_local = "numeric",
                            individuals_meta = "numeric",
                            species_meta = "numeric",
                            speciation_local = "numeric",
                            speciation_meta = "numeric",
                            extinction_meta = "numeric",
                            trait_sigma = "numeric",
                            env_sigma = "numeric",
                            comp_sigma = "numeric",
                            dispersal_prob = "numeric",
                            mutation_rate = "numeric" ,
                            equilib_escape = "numeric",
                            num_basepairs = "numeric",
                            init_type = "character"
                          ))

roleParams <- function(){
  return(new("roleParams"))
}