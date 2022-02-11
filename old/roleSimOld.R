roleSimOld <- function(params = NULL, startModel = NULL, initType = "oceanic_island", niter = 100, nruns = 1, niter_timestep = 50, print = FALSE) {
  
  # if params SPECIFIED and startModel UNSPECIFIED, get niter and nruns from params
  if((!is.null(params) & is.null(startModel)) | (!is.null(params) & is.null(startModel)))
  {
    nruns <- params@nruns
    niter <- params@niter
  }
  # else params is UNSPECIFIED and startModel SPECIFIED get params from startModel
  else if(is.null(params) & is.null(startModel))
  { 
    params <- startModel@params
    params <- roleParams(nruns,niter,niter_timestep)
    setDefaultParams(params)
  }
  # if params and startModel specified, use params
  if(!is.null(params) && !is.null(startModel))
  {
    
  }
  
  # convert params and init to C++ objects
  #params <- toCpp(params)
  
  # NOTE - currently if starting with an existing sim, the params of that sim are used and not the new params
  if(is.null(startModel)) { # if no init sim to start with, initialize a new sim with params
    init <- initSim(params)
    init$print <- print 
  }
  else{
    # convert model to C++ model
    init <- toCpp(startModel)
  }
  
  # initialize list to hold finished model runs
  out <- list()
  # do nrun model runs 
  for(i in 1:nsim)
  {
    # set params for new run
    init$params <- params@values[[1]]
    
    # iterate over the simulation and return the output
    out[[i]] <- iterSim(init, niter, niter_timestep, print) # else add the next sim to the list 
  }
  
  # return list of sims
  return(out)
}

# ----
#' @description function to initialize role simulation, returns a `roleModelCpp`
#' object
#' @param params a roleParamsCpp object containing the model parameters - if null, defaults are set

initSim <- function(params, type) {
  
  # simulate phylogeny
  phy <- TreeSim::sim.bd.taxa(params@values[["species_meta"]][[1]], numbsim = 1,
                              lambda = params@values[["speciation_meta"]][[1]],
                              mu = params@values[["extinction_meta"]][[1]], complete = FALSE)[[1]]
  
  # simulate metacommunity SAD
  abundance_m <- .lseriesFromSN(params@values[["species_meta"]][[1]],
                                params@values[["individuals_meta"]][[1]])
  
  # initalize meta comm traits:
  # first column is species ID, second column is trait value
  # Smax rows
  traits_m <- cbind(1:params@values[["species_meta"]][[1]],
                    ape::rTraitCont(phy, sigma = params@values[["trait_sigma"]][[1]]))
  
  # create metaCommCpp object
  meta <- new(metaCommCpp, abundance_m, traits_m, params@values[["species_meta"]][[1]])
  
  # vector of 0 abundances
  abundance_l <- rep(0, params$values$species_meta * 100)
  
  # initialize local species abundances with one species having all individuals
  # eventually may want different ways of initializing abundances based on Harmon paper
  # argument in params  specifying abundance initalization model
  # save sampling index and use to assign to correct species of local abundances,
  # and the appropriate trait from the metacomm pool and add to trait matrix
  
  if(type == "oceanic_island"){
    # index of the species that will initially have all abundance
    i <- sample(params$values$species_meta, 1, prob = meta$abundance)
    # passing all abundance to that species
    abundance_l[i] <- params$values$individuals_local
  }
  
  else if(type == "bridge_island"){
    abundance_l <- sample(params$values$species_meta, prob = meta$abundance)
  }
  
  # counter keeping track of max number of possible species in local comm
  Smax_ <- params$values$species_meta
  
  # init local species traits
  traits_l <- numeric(100)
  
  # extract trait for starting species from meta
  traits_l[i] <- meta$traits[which(meta$traits[meta$traits[,1]] == i),2]
  
  # soon set pi to be simulated sequences, and update outside of C++ over time
  pi_l <- rep(1:Smax_)
  
  # create localCommCpp object
  local <- new(localCommCpp, abundance_l, traits_l, Smax_, pi_l)
  
  # convert ape phylo to rolePhylo
  phy <- apeToRolePhylo(phy)
  
  # convert rolePhylo to rolePhyloCpp
  phy <- rolePhyloToCpp(phy)
  
  # create roleModelCpp object of local comm, meta comm, phylogeny, and params
  out <- new(roleModelCpp,local,meta,phy,params)
  
  return(out)
}