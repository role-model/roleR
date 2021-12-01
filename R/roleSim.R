#' @title RoLE model simulation
#'
#' @description Simulate communities under the RoLE model. The key distinction
#' between the two functions is that \code{roleSim} is optimized to run many
#' simulations, while \code{roleSimPlay} is meant to run one simulation with
#' periodic output of results for visualization and exploration.
#' \code{roleSimPlay} is intended primarily for the accompanying Shiny App. R
#'
#'
#' @param params list of parameters
#' @param init an optional initial condition specified as a \code{roleModelCpp}
#' object
#' @param nstep number of simulation steps to run; if \code{prop_equilib} is
#' specified in the \code{params} list, \code{nstep} will over-ride that
#' parameter
#' @param nsim number of simulations to run
#' @param nout frequency of intermediate results output
#'
#' @details Stub
#'
#' @rdname roleSim
#' @export

roleSim <- function(params, init = NULL, nstep = NULL, nsim = 1) {

    # initialize simulation
    if(is.null(init)) {
        init <- initSim(params)
    }

    # loop over steps
    iterSimCpp(init, nstep)
}

# ----
#' @description function to initialize role simulation, returns a `roleModelCpp`
#' object
#' @param params a roleParamsCpp object containing the model parameters - if null, defaults are set

initSim <- function(params = NULL) {

    # if no roleParamsCpp object provided
    if(is.null(params))
    {
        #create empty paramValuesCpp object
        #empty constructor includes a plausible default set of values
        vals <- new(paramValuesCpp)

        #create roleParamsCpp object from values
        params <- new(roleParamsCpp,vals,"sim", 1)
    }

    # simulate phylogeny
    phy <- TreeSim::sim.bd.taxa(params$values$species_meta, numbsim = 1,
                                lambda = params$values$speciation_meta,
                                mu = params$values$extinction_meta, complete = FALSE)[[1]]

    # simulate metacommunity SAD
    abundance_m <- .lseriesFromSN(params$values$species_meta,
                                params$values$individuals_meta)

    # initalize meta comm traits:
    # first column is species ID, second column is trait value
    # Smax rows
    traits_m <- cbind(1:params$values$species_meta,
                         ape::rTraitCont(phy, sigma = params$values$trait_sigma))

    # create metaCommCpp object
    meta <- new(metaCommCpp, abundance_m, traits_m, params$values$species_meta)

    # vector of 0 abundances
    abundance_l <- rep(0, params$values$species_meta * 100)

    # initialize local species abundances with one species having all individuals
    # eventually may want different ways of initializing abundances based on Harmon paper
    # argument in params  specifying abundance initalization model
    # save sampling index and use to assign to correct species of local abundances,
    # and the appropriate trait from the metacomm pool and add to trait matrix

    # index of the species that will initially have all abundance
    i <- sample(params$values$species_meta, 1, prob = meta$abundance)

    # passing all abundance to that species
    abundance_l[i] <- params$values$individuals_local

    # counter keeping track of max number of possible species in local comm
    Smax_ <- params$values$species_meta

    # extract local comm traits from meta traits
    traits_l <- matrix(NA, nrow = params$values$species_meta * 10000, ncol = 2)

    #changed this - correct?
    traits_l[i,] <- c(i, meta$traits[meta$traits[,1] == i, 2])

    pi_l <- rep(1:Smax_)

    # create localCommCpp object
    local <- new(localCommCpp, abundance_l, traits_l, Smax_,pi_l)

    # convert ape phylo to rolePhylo
    #phy <- as(phy, "rolePhylo")
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

.lseriesFromSN <- function(S, N) {
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

rolePhyloToCpp <- function(phylo){
    n <- phylo@n

    e <- phylo@e
    l <- phylo@l
    alive <- phylo@alive
    tipNames <- phylo@tipNames
    scale <- phylo@scale
    out <- new(rolePhyloCpp,n,e,l,alive,tipNames,scale)
    return(out)
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
