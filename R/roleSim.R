#' @title An S4 class containing multiple model runs and the params used to generate them
#'
#' @slot modelRuns a list of roleModel objects, each containing its own current state and 
#' a timeseries of past states 
#' @slot params a roleParams object
#' @slot rootSim 
#'
#' @export

setClass('roleSim',
         slots = c(modelRuns = 'list',
                   params = 'roleParams'))

# NOTE - should we require users to specify a params object, and thus avoid needing niter, nruns, niter_timestep as params of roleSim? 
#' @title Run a RoLE model simulation and return a roleSim object
#'
#' @description Simulate communities under the RoLE model, The key distinction
#' between the two functions is that \code{roleSim} is optimized to run many
#' simulations, while \code{roleSimPlay} is meant to run one simulation with
#' periodic output of results for visualization and exploration.
#' \code{roleSimPlay} is intended primarily for the accompanying R Shiny App
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
#' sim <- roleSim(params)
#' 
#' sim <- roleSim(params, startModel = sim@runs[1])
#' 
#' @return a roleSim object
#'
#' @rdname roleSim
#' @export

roleSim <- function(params, startModel = NULL, initType = "oceanic_island", print = FALSE) {
  
    nruns <- params@nruns
    niter <- params@niter
    niter_timestep <- params@niter_timestep
    
    # prepare params for C++ loop by stretching 1-length values to niter and sampling from provided priors 
    cparams <- stretchAndSampleParams(params)
    
    # initialize a C++ roleModelCpp starting condition
    init <- initSim(cparams,initType) 
    
    # initialize list to hold finished model runs
    runs <- list()
    
    # do nrun model runs 
    for(i in 1:nruns)
    {
      # set the named list of param values
      init@params <- cparams@values[[1]]
      
      # iterate over the simulation and return the output
      runs[[i]] <- roleModelFromCpp(iterSim(init, niter, niter_timestep, print)) # else add the next sim to the list 
    }
    
    # return list of sims
    return(new("roleSim", runs, params))
}

roleSimPlay <- function() {
}

# ----
#' @description function to initialize role simulation, returns a `roleModelCpp`
#' object
#' @param params a roleParamsCpp object containing the model parameters - if null, defaults are set

initSim <- function(params, type) {
    
    # some testing lines
    #params <- roleParams(nrun=1,niter=100,niter_timestep=10,defaults=TRUE)
    #params <- setDefaultParams(params)
    
    # simulate phylogeny
    phy <- TreeSim::sim.bd.taxa(params@values[[1]][["species_meta"]], numbsim = 1,
                                lambda = params@values[[1]][["speciation_meta"]],
                                mu = params@values[[1]][["extinction_meta"]], complete = FALSE)[[1]]

    # simulate metacommunity SAD
    abundance_m <- lseriesFromSN(params@values[[1]][["species_meta"]],
                                  params@values[[1]][["individuals_meta"]])
    
    # initalize meta comm traits:
    # first column is species ID, second column is trait value
    # Smax rows
    traits_m <- cbind(1:params@values[[1]][["species_meta"]],
                         ape::rTraitCont(phy, sigma = params@values[[1]][["trait_sigma"]]))

    # create metaCommCpp object
    meta <- new(metaCommCpp, abundance_m, traits_m, params@values[[1]][["species_meta"]])
  
    
    # initialize vector of 0 species abundances
    abundance_l_sp <- rep(0, params@values[[1]][["species_meta"]])
    
    # oceanic island model assigns all abundance to one species
    if(type == "oceanic_island"){
      # index of the species that will initially have all abundance
      i <- sample(params@values[[1]][["species_meta"]], 1, prob = meta$abundance)
      # passing all abundance to that species
      abundance_l_sp[i] <- params@values[[1]][["individuals_local"]]
    }
    # bridge island model assigns abundances to all species proportional to species abundance
    else if(type == "bridge_island"){
      # vector of species
      abundance_l_sp <- sample(params@values[[1]][["species_meta"]], prob = meta$abundance)
    }
    
    # counter keeping track of max number of possible species in local comm
    Smax_ <- params@values[[1]][["species_meta"]]

    # init local species traits
    traits_l_sp <- traits_m
    
    # soon set pi to be simulated sequences
    pi_l <- rep(1:Smax_)
    
    # create localCommCpp object
    local <- new(localCommCpp, abundance_l_sp, traits_l_sp, Smax_, pi_l)

    # convert ape phylo to rolePhylo
    phy <- apeToRolePhylo(phy)

    # convert rolePhylo to rolePhyloCpp
    phy <- rolePhyloToCpp(phy)

    # create roleModelCpp object of local comm, meta comm, phylogeny, and params
    out <- new(roleModelCpp,local,meta,phy,params)

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
#' @param model object of class \code{roleSim}
#' @param dir the directory to write to i.e. "data/models" 
#' @param fileName the name of the new file
#' @param saveTxt specifies whether to save a .txt describing some sim information
#' for viewing outside of R}
#'
#' @export

writeModel <- function(model, dir = NULL, fileName, saveTxt = TRUE)
{
  if(is.null(dir)){
    dir = getwd()
  }
  saveRDS(model,paste0(dir,"/",fileName,".role"))
  if(saveTxt){
    con <- file(paste0(dir, "/", fileName, "_info", ".txt"))
    writeLines(c("Info for ROLE model object","Created by Jacob Idec"), con)
    close(con)
  }
}
