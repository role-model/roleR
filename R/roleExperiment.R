#' @title An S4 class containing multiple model runs and the params used to generate them
#'
#' @slot modelRuns a list of roleModel objects, each containing its own current state and 
#' a timeseries of past states 
#' @slot params a roleParams object
#'
#' @export

setClass('roleExperiment',
         slots = c(modelRuns = 'list',
                   params = 'roleParams'))

# NOTE - should we require users to specify a params object, and thus avoid needing niter, nruns, niter_timestep as params of roleExperiment? 
#' @title Run a RoLE model simulation and return a roleExperiment object
#'
#' @description Simulate communities under the RoLE model, The key distinction
#' between the two functions is that \code{roleExperiment} is optimized to run many
#' simulations, while \code{roleExperimentPlay} is meant to run one simulation with
#' periodic output of results for visualization and exploration.
#' \code{roleExperimentPlay} is intended primarily for the accompanying R Shiny App
#'
#' @param params a \code{roleParams} object - if left unspecified a default set of "plausible"
#' parameter values is used 
#' @param startModel a roleModel object - if unspecified is built from scratch using params and initType
#' @param initType a string specifying the type of model to create if starting from scratch without a startModel -
#' either "oceanic_island" or "bridge_island"
#' @param print a bool specifying whether to print some information as the model runs
#' 
#' @example 
#' params <- roleParams(niter = 5000,defaults=TRUE) 
#' role <- roleExperiment(params)
#' 
#' role <- roleExperiment(params, startModel = role@runs[1])
#' 
#' @return a roleExperiment object
#'
#' @rdname roleExperiment
#' @export

roleExperimentPlay <- function()
{
  for(t in timesteps)
  {
    iterSim()
    # write sim 
  }
}


roleExperiment <- function(params, startModel = NULL, initType = "oceanic_island", print = FALSE) {
  
    nruns <- params@nruns
    niter <- params@niter
    niter_timestep <- params@niterTimestep
    
    # prepare params for C++ loop by stretching 1-length values to niter and sampling from provided priors 
    cparams <- stretchAndSampleParams(params)
    
    # initialize list to hold finished model runs
    runs <- list()
    
    # do nrun model runs 
    for(i in 1:nruns)
    {
      i <- 1
      # get the params for the run 
      parlist <- cparams@values[[i]] 
      
      # initialize a C++ roleModelCpp starting condition
      init <- initModel(parlist,initType, niter) 
      
      init$print <- print
      
      # iterate over the simulation and return the output
      runs[[i]] <- roleModelFromCpp(iterSim(init, niter, niter_timestep, print)) # else add the next sim to the list 
    }
    
    # return list of sims
    return(new("roleExperiment", runs, params))
}

roleExperimentPlay <- function() {
}

# ----
#' @description function to initialize role simulation, returns a `roleModelCpp`
#' object
#' @param params a named List of parameters 

initModel <- function(parlist, type, niter) {
    
    # set aug length
    naug <- niter  
  
    # some testing lines
    params <- roleParams(nrun=1,niter=100,niterTimestep=10,defaults=TRUE)

    # simulate phylogeny
    phy <- TreeSim::sim.bd.taxa(parlist[["species_meta"]][1], numbsim = 1,
                                lambda = parlist[["speciation_meta"]][1],
                                mu = parlist[["extinction_meta"]][1], complete = FALSE)[[1]]
    
    # simulate metacommunity SAD
    abundance_m <- lseriesFromSN(parlist[["species_meta"]][1],
                                  parlist[["individuals_meta"]][1])
    
    # initalize meta comm traits:
    # first column is species ID, second column is trait value
    # Smax rows
    # talk over with andy - why is rTraitCont producing values in the thousands and negative thousands? 
    traits_m <- cbind(1:parlist[["species_meta"]][1],
                         ape::rTraitCont(phy, sigma = parlist[["trait_sigma"]][1]))
    # create metaCommCpp object
    meta <- new(metaCommCpp, abundance_m, traits_m)
    
    # initialize vector of 0 species abundances
    abundance_l_sp <- rep(0, parlist[["species_meta"]][1])
    
    # oceanic island model assigns all abundance to one species
    if(type == "oceanic_island"){
      # index of the species that will initially have all abundance
      i <- sample(parlist[["species_meta"]][1], 1, prob = meta$abundance)
      # passing all abundance to that species
      abundance_l_sp[i] <- parlist[["individuals_local"]][1]
    }
    # bridge island model assigns abundances to all species proportional to species abundance
    else if(type == "bridge_island"){
      # vector of species
      abundance_l_sp <- sample(parlist[["species_meta"]][1], prob = meta$abundance)
    }
    

    # init local species traits
    traits_l_sp <- traits_m
    
    # soon set pi to be simulated sequences (of length equal to the number of species in the metacomm)
    pi_l <- rep(1:parlist[["species_meta"]][1])
    
    # create localCommCpp object
    local <- new(localCommCpp, abundance_l_sp, traits_l_sp, pi_l, naug)

    # convert ape phylo to rolePhylo
    phy <- apeToRolePhylo(phy)
  
    # convert rolePhylo to rolePhyloCpp
    phy <- rolePhyloToCpp(phy)
    
    # create roleModelCpp object of local comm, meta comm, phylogeny, and params
    out <- new(roleModelCpp,local,meta,phy,parlist)

    return(out)
}

# ----
#' @description function to solve for parameter of logseries
#' @param S number of species
#' @param N number of individuals

lseriesFromSN <- function(S, N) {
    # solve for alpha paramter
    asol <- uniroot(interval = c(.Machine$double.eps^0.25,
                                 .Machine$integer.max),
                    f = function(a) {
                        a * log(1 + N / a) - S
                    })

    # calculate p parameter and beta (as used by pika)
    p <- 1 - exp(-S / asol$root)
    beta <- -log(p)

    # calculate idealized SAD from parameter
    thisSAD <- pika::sad(model = 'fish', par = beta)
    thisSAD <- pika::sad2Rank(thisSAD, S = S)

    # return relative abundances
    return(thisSAD / sum(thisSAD))
}

apeToRolePhylo <- function(phylo){
  
  # extract number of times
  n <- ape::Ntip(phylo)
  
  # extract edge matrix and edge lengths
  e <- phylo$edge
  l <- phylo$edge.length
  
  # extract tip labels
  tipNames <- phylo$tip.label
  
  # calculate alive or not
  tipAge <- ape::node.depth.edgelength(phylo)[1:n]
  
  alive <- rep(TRUE, n)
  alive[tipAge < max(tipAge)] <- FALSE
  
  # set default scale
  scale <- 1
  
  return(rolePhylo(n = n, e = e, l = l, alive = alive,
                   tipNames = tipNames, scale = scale))
}

apeToRolePhyloOld <- function(phylo){

    # extract number of times
    n <- ape::Ntip(phylo)

    # extract edge matrix and edge lengths
    e <- phylo$edge
    l <- phylo$edge.length

    # extract tip labels
    tipNames <- phylo$tip.label

    # calculate alive or not
    tipAge <- ape::node.depth.edgelength(phylo)[1:n]

    alive <- rep(TRUE, n)
    alive[tipAge < max(tipAge)] <- FALSE

    # set default scale
    scale <- 1

    # buffer objects so we can add new species without augmenting objects
    addOn <- n * 100
    e <- rbind(e, matrix(-1, nrow = addOn, ncol = 2))
    l <- c(l, rep(0, addOn))
    alive <- c(alive, rep(FALSE, addOn))
    tipNames <- c(tipNames, rep('', addOn))

    return(rolePhylo(n = n, e = e, l = l, alive = alive,
                     tipNames = tipNames, scale = scale))
}

apeToPhyloCpp <- function(phylo){
    n <- ape::Ntip(phylo)

    e <- phylo$edge
    l <- phylo$edge.length
    tipNames <- phylo$tip.label
    tipAge <- ape::node.depth.edgelength(phylo)[1:n]
    alive <- rep(TRUE,n)
    alive[tipAge < max(tipAge)] <- FALSE;
    scale <- 1;
    out <- new(rolePhyloCpp,n,e,l,alive,tipNames,scale)
    return(out)
}

#' @title Write a ROLE sim object to a serialized R file for transfer, sharing, 
#' or to use with Docker. File contents can only be read within R. 
#'
#' @param model object of class \code{roleExperiment}
#' @param dir the directory to write to i.e. "data/models" 
#' @param fileName the name of the new file
#' @param saveTxt specifies whether to save a .txt describing some sim information
#' for viewing outside of R}
#'
#' @export

writeRoleExperiment <- function(role, dir = NULL, fileName, saveTxt = TRUE)
{
  if(is.null(dir)){
    dir = getwd()
  }
  saveRDS(role,paste0(dir,"/",fileName,".roleexperiment"))
  if(saveTxt){
    con <- file(paste0(dir, "/", fileName, ".roleexperimentinfo", ".txt"))
    title_line = "Metadata for RoLE experiment R object - load into R using readRDS(file_location_and_name)"
    info_line = paste("Contains", role@params@nruns, "runs of", role@params@niter, "iterations")
    author_line = paste("Author:", role@params@meta[1])
    date_line = paste("Date:", role@params@meta[2])
    desc_line = paste("Description:", role@params@meta[3])
    writeLines(c(title_line,info_line,author_line,date_line,desc_line), con)
    close(con)
  }
}

dummyExperiment <- function(run=FALSE)
{
  params <- roleParams(nrun=1,niter=1000,niterTimestep=100,defaults=TRUE)
  cparams <- stretchAndSampleParams(params)
  parlist <- cparams@values[[1]]
  model <- initModel(parlist,type="bridge_island",niter=1000)
  
  model$print <- FALSE
  model$local$print <- FALSE
  model$timeseries[[1]] <- model$copyData(1)
  
  if(run){
    iterSim(model,params@niter,params@niterTimestep,FALSE)
  }
  model <- roleModelFromCpp(model)
  runs <- list(model) 
  
  # return list of sims
  return(new("roleExperiment", modelRuns=runs, params=params))
}