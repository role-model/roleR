#' @title An S4 class containing parameters for a roleExperiment
#' 
#' @details Values is an array of named lists, 1 for each roleModel run of roleExperiment. 
#' Each named list contains 12 numeric vectors (1 per parameter) with each vector being either a single value to be held constant across all
#' iterations of the model run OR an niter length vector specifying the value of that parameter for each individual iteration
#' of the model run 
#'
#' @slot nruns an integer specifying the number of runs (number of roleModels) to simulate -  defaults to 1 
#' @slot niter an integer specifying the number of iterations per model run 
#' @slot niterTimestep an integer specifying the frequency at which step states for each model are saved in the model timeseries 
#' i.e. a value of 5 means save every 5 iteration steps
#' @slot values an array of named lists of parameter values
#' @slot priors priors to draw values from when the model starts - any priors specified OVERRIDE their params contained in values
#' @slot iterfuncs functions to draw values from which take an iteration step and return a value
#' 
#' @export

roleParams <- setClass('roleParams',
         slots = c(
           nruns = "numeric", 
           niter = "numeric",
           niterTimestep = "numeric",
           values = "array",
           priors = "array",
           iterfuncs = "array",
           meta = "character"
         ), 
         
         validity=function(object)
         {
           # make sure nruns and niter are valid numbers
           if(object@nruns <= 0){
             return(FALSE)
           }
           if(object@niter <= 1){
             return(FALSE)
           }
           
           # check all params to make sure they are either 1 length or niter length
           # for each run
           for(r in 1:length(object@values))
           {
             # for each param
             for(p in 1:length(object@values[[1]]))
             {
               # return false if not 1-length or niter length
               if((length(object@values[[r]][[p]]) != 1) & (length(object@values[[r]][[p]]) != object@niter)){
                 # AND there are no priors for that param name
                 if(length(object@priors[[1]][[p]]) != 1){
                  return(FALSE)}
               }
             }
           }
           
           return(TRUE)
         })

#' @title Create a new roleParams object
#'
#' @param nruns an integer specifying the number of runs (number of roleModels) to simulate -  defaults to 1 
#' @param niter an integer specifying the number of iterations per model run 
#' @param niterTimestep an integer specifying the frequency at which step states for each model are saved in the model timeseries 
#' i.e. a value of 5 means save every 5 iteration steps
#' @param defaults a bool specifying whether to initially assign default values to params 
#' @param individuals_local a value or set of values to assign to params
#'
#' @examples 
#' params <- roleParams(nruns = 3,niter = 10000, niter_timestep = 25)
#' 
#' @export

roleParams <- function(nruns, niter, niterTimestep=10,defaults=TRUE,
                       
                       individuals_local=NA,individuals_meta=NA, 
                       species_meta=NA,speciation_local=NA,
                       speciation_meta=NA, extinction_meta=NA,
                       trait_sigma=NA,env_sigma = NA,
                       comp_sigma=NA,dispersal_prob=NA,
                       mutation_rate = NA, equilib_escape=NA,
                       genesim_savesteps=NA,num_basepairs=NA){

  values <- array(list(),nruns)
  
  # set default zeroes
  # could also set default values but I dunno
  for(i in 1:nruns)
  {
    values[[i]][["individuals_local"]] <- individuals_local
    values[[i]][["individuals_meta"]] <- individuals_meta
    values[[i]][["species_meta"]] <- species_meta
    values[[i]][["speciation_local"]] <- speciation_local
    values[[i]][["speciation_meta"]] <- speciation_meta
    values[[i]][["extinction_meta"]] <- extinction_meta
    values[[i]][["trait_sigma"]] <- trait_sigma
    values[[i]][["env_sigma"]] <- env_sigma
    values[[i]][["comp_sigma"]] <- comp_sigma
    values[[i]][["dispersal_prob"]] <- dispersal_prob
    values[[i]][["mutation_rate"]] <- mutation_rate
    values[[i]][["equilib_escape"]] <- equilib_escape
    values[[i]][["genesim_savesteps"]] <- genesim_savesteps
    values[[i]][["num_basepairs"]] <- num_basepairs
  }
  
  if(is.null(niterTimestep)){
    niterTimestep = 10}

  params <- new("roleParams", nruns=nruns, niter=niter,values=values,niterTimestep=niterTimestep)
  if(defaults){
    params <- setDefaultParams(params)}
  
  return(params)
}


#' @title Set a param of roleParams
#'
#' @param x a roleParams object
#' @param paramName a string specifying the name of the parameter to set 
#' @param values a vector containing one parameter value or a timeseries of parameter values
#' @param runNum an integer specifying 
#' 
#' @examples 
#' setParam(params, paramName = "speciation_local", values = 0.125)
#' setParam(params, paramName = "mutation_rate", values = 0.2, runNum = 3)
#' vals <- my_value_generating_function(niter) # NOTE - is this a good use case? are there cases where people won't want to set a prior function and instead create a series and feed it in? I think if using a series of complex functions then yes 
#' setParam(params, paramName = "mutation_rate", values = vals)
#' 
#' NOTE -  may want to permit setting values starting at a specified iter n 
#' @export

setGeneric('setParam', function(x, ...) standardGeneric('setParam'), signature = 'x')
setMethod("setParam", signature(x="roleParams"),
          function(x,paramName,values,runNum=NULL) {
            
            # if no runNum provided, apply to all sims
            if(is.null(runNum))
            {
              for(i in 1:x@nruns)
              {
                x@values[[i]][[paramName]] <- values
              }
            }
            
            # else apply only to the target sim_num
            else
            {
                x@values[[runNum]][[paramName]] <- values
            }
            
            return(x)
          }
)

#' @title Set a prior function to use when sampling parameters of roleParams
#'
#' @details if a prior is specified, it overrides any previous param values

#' @param x a roleParams object
#' @param paramName a string specifying the name of the parameter to set a prior for
#' @param priorFunction an R function to use when sampling i.e. runif(min=0,max=1) 
#' @param runNum an integer specifying the runNum to set a prior for
#' NOTE -  may want to permit setting values starting at a specified iter n 
#' 
#' @examples 
#' setPrior(params, paramName = "speciation_local", priorFunction = runif(min=0,max=0.25))
#' 
#' @export

setGeneric('setPrior', function(x, ...) standardGeneric('setPrior'), signature = 'x')
setMethod("setPrior", signature(x="roleParams"),
          function(x,paramName,priorFunction,runNum=NULL) {
            
            # if no runNum provided, apply to all sims
            if(is.null(runNum))
            {
              for(i in 1:x@nruns)
              {
                x@priors[[i]][[param_name]] <- priorFunction
              }
            }
            
            # else apply only to the target sim_num
            else
            {
              x@priors[[runNum]][[param_name]] <- priorFunction
            }
            
          return(x)
          }
)



#' @title Set a default set of params in roleParams
#'
#' @param x a roleParams object
#' 
#' @examples 
#' setDefaultParams(params)
#' 
#' @export

setGeneric('setDefaultParams', function(x, ...) standardGeneric('setDefaultParams'), signature = 'x')
setMethod("setDefaultParams", signature(x="roleParams"),
          function(x) {
            x <- setParam(x,"species_meta",15)
            x <- setParam(x,"individuals_meta",500)
            x <- setParam(x,"individuals_local",100)
            x <- setParam(x,"dispersal_prob",0.5)
            x <- setParam(x,"speciation_local",0.1)
            x <- setParam(x,"extinction_meta",0.1)
            x <- setParam(x,"speciation_meta",0.1)
            x <- setParam(x,"trait_sigma",0.5)
            x <- setParam(x,"env_sigma",0.1)
            x <- setParam(x,"comp_sigma",0.1)
            x <- setParam(x,"trait_z",0.5)
            x <- setParam(x,"genesim_nstep",10)
            return(x)
          }
)

# prep the params to be read into C++ - use prior distributions if they are set, and stretch 1 value params across all iters
# NOTE -  only used in roleExperiment method so will ultimately be dotted and inaccessible to user
stretchAndSampleParams <- function(params) {
  
  # if priors specified, sample params using priors
  # if 1 set of priors specified
  if(length(params@priors) == 1)
  {
    runPriors <- params@priors[[1]]
    for(j in 1:length(runPriors))
    {
      # get function
      funName <- names(runPriors)[j]
      fun <- runPriors[[funName]]
      # get values using fun with n = niter, and assign to out params
      params@values[[i]][[param_name]] <- fun(n=params@niter)
    }
  }
  
  # TODO - if multiple sets of priors specified for each model run  
  else if(length(params@priors) > 1){
    for(i in 1:length(params@priors))
    {
    }
  }
  
  # stretch 1-value params across all iter values
  # for each run
  for(r in 1:length(params@values))
  {
    # for each param
    for(p in 1:length(params@values[[1]]))
    {
      if(length(params@values[[r]][[p]]) == 1){
        params@values[[r]][[p]] <- rep(params@values[[r]][[p]],params@niter)
      }
    }
  }
  
  return(params)
}

writeRoleParams <- function(params, dir = NULL, fileName, saveTxt = TRUE)
{
  if(is.null(dir)){
    dir = getwd()
  }
  saveRDS(params,paste0(dir,"/",fileName,".roleparams"))
  if(saveTxt){
    con <- file(paste0(dir, "/", fileName, ".roleparamsinfo", ".txt"))
    title_line = "Info: helper metadata for RoLE params R object - load object into R using readRDS(file_location_and_name)"
    info_line = paste("Contains", params@nruns, "runs of", params@niter, "iterations")
    author_line = paste("Author:", params@meta[1])
    date_line = paste("Date:", params@meta[2])
    desc_line = paste("Description:", params@meta[3])
    writeLines(c(title_line,info_line,author_line,date_line,desc_line), con)
    close(con)
  }
}