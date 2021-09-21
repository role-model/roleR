#' @title RoLE model simulation
#'
#' @description Simulate communities under the RoLE model. The key distinction
#' between the two functions is that \code{roleSim} is optimized to run many
#' simulations, while \code{roleSimPlay} is meant to run one simulation with
#' periodic output of results for visualization and exploration.
#' \code{roleSimPlay} is intended primarily for the accompanying Shiny App.
#'
#'
#' @param params list of parameters
#' @param init an optional initial condition specified as a \code{roleComm}
#' object
#' @param nstep number of simulation steps to run; if \code{prop_equilib} is
#' specified in the \code{params} list, \code{nstep} will over-ride that
#' parameter
#' @param nsim number of simulations to run
#' @param nout frequency of intermediate results output
#'
#' @details Stub
#'
#' @examples
#' params <- list(species_meta = 100,
#'                individuals_meta = 10000,
#'                individuals_local = 1000,
#'                dispersal_prob = 0.1,
#'                speciation_local = 0.01)
#' testSim <- roleSim(params, nstep = 10, nsim = 1)
#' testSim
#'
#' @return An object of class \code{roleComm} with the following elements:
#' \describe{
#'   \item{\code{local_comm}}{a list with three named elements:
#'   \code{Abundance}, \code{Traits}, \code{pi}}
#'   \item{\code{meta_comm}}{a list with two named elements: \code{Abundance},
#'   \code{Traits}}
#'   \item{\code{phylo}}{a phylogeny (of class \code{ape::phylo}) for all
#'   species (extinct and extant) in the meta and local communities}
#'   \item{\code{params}}{a named list of model parameters}
#' }
#'
#' In the case of \code{roleSimPlay}, a list of \code{roleComm} objects is
#' returned
#'
#' @rdname roleSim
#' @export

roleSim <- function(params, init = NULL, nstep = NULL, nsim = 1) {
    # ----
    # initialize simulation
    if(is.null(init)) {
        init <- .initSim(params)
    }

    # make sure `params` reflect the user-specified values (potentially
    # overwriting old values)
    init@params <- params

    # ----
    # loop over steps
    .iterSim(init, nstep)
}

#' @rdname roleSim
#' @export
#'
roleSimPlay <- function(params, init = NULL, nstep = NULL, nout = NULL) {
    # ----
    # initialize simulation

    if(is.null(init)) {
        init <- .initSim(params)
    }

    # make sure `params` reflect the user-specified values (potentially
    # overwriting old values)
    init@params <- params

    # ----
    # loop over steps

    # number of chunks to iterate over
    B <- floor(nstep / nout)
    if(B <= 1) B <- 2
    allSims <- vector('list', B)

    # first set of iterations
    allSims[[1]] <- .iterSim(init, nout)

    # loop over remaining iterations
    for(b in 2:B) {
        # plotting stuff, etc goes here

        # update `roleComm` object
        allSims[[b]] <- .iterSim(allSims[[b - 1]], nout)
    }

    # return all the simulations
    return(allSims)
}

# ----
#' @description function to initialize role simulation, returns a `roleComm`
#' object
#' @param params list of model parameters

.initSim <- function(params) {

    p <- params

    #if null params, initialize default plausible set of params
    # change to set any non input to defaults
    if(is.null(p))
        p <- list(species_meta = 15,
                   individuals_meta = 500,
                   individuals_local = 100,
                   dispersal_prob = 0.5,
                   speciation_local = 0.1,
                   extinction_meta = 0.8,
                   speciation_meta = 1,
                   trait_sigma = 0.1)

    #create roleParams object
    p <- roleParams(p,"sim")

    # simulate phylogeny
    phy <- TreeSim::sim.bd.taxa(p@params$species_meta, numbsim = 1,
                                lambda = p@params$speciation_meta,
                                mu = p@params$extinction_meta)[[1]]

    # create metaComm object
    meta <- metaComm(abundance = numeric(), traits = matrix(numeric()),
                     Smax = 0)
    # simulate metacommunity SAD
    meta@abundance <- .lseriesFromSN(p@params$species_meta,
                                             p@params$individuals_meta)
    meta@Smax <- p@params$species_meta

    #initalize meta comm traits without rep NA
    local@traits <- matrix(c(ape::rTraitCont(phy, sigma = p@params$trait_sigma),
                             rep(NA, p@params$species_meta * 99)))

    # create localComm object
    local <- localComm(abundance = numeric(),
                       traits = matrix(numeric()), pi = numeric(), Smax = 0)
    # vector of 0 abundances
    local@abundance <- rep(0, p@params$species_meta * 100)

    # initialize local species abundances with one species having all individuals
    # eventually may want different ways of initializing abundances based on Harmon paper
    # argument in params  specifying abundance initalization model
    # save sampling index and use to assign to correct species of local abundances,
    # and the appropriate trait from the metacomm pool and add to trait matrix

    local@abundance[sample(p@params$species_meta, 1,
                                   prob = meta@abundance)] <-
        p@params$individuals_local

    # counter keeping track of max number of possible species in local comm
    local@Smax <- p@params$species_meta

    # simulate local comm traits
    local@traits <- matrix(c(ape::rTraitCont(phy, sigma = p@params$trait_sigma),
                              rep(NA, p@params$species_meta * 99)))

    # convert ape phylo to rolePhylo
    phy <- as(phy, "rolePhylo")

    # create roleModel object of local comm, meta comm, phylogeny, and params
    out <- roleModel(local,meta,phy,p)

    # General note on book keeping: for vectors that we expect to grow (i.e.,
    # traits and local abundances) we make a vector 100 times longer than the
    # initial number of species in the metacommunity rather than needing to
    # augment those vectors

    return(out)
}

# Altered to fit new methods but incompatible with old initSim
# ----
# function to iterate role simulation
#' @param model an object of class `roleModel`
#' @param nstep number of steps to iterate

.iterSim <- function(model, nstep) {
    # browser()
    for(i in 1:nstep) {

        # calls deathRole, which samples a species for death then calls death on localComm & rolePhylo
        model <- death(model)

        # if speciation occurs based on speciation chance param
        # local speciation param may be renamed
        if(runif(1) <= model@params@speciation_local) {

            # call speciation on localComm & rolePhylo
            model <- speciation(model)
            # model@localComm <- speciation(model@localComm)
            # model@rolePhylo <- speciation(model@rolePhylo)
        }

        # if immigration occurs based on dispersal chance param
        else if(runif(1) <= model@params@dispersal_prob) {

            # sample a species for immigration relative to metacomm abundance
            i <- sample(model@params@species_meta, 1,
                          prob = model@metaComm@abundance)

            # call immigration on localComm
            model@localComm <- immigration(model@localComm)

        }

        # else birth occurs
        else {
            # sample a species for birth relative to local abundance
            i <- sample(model@localComm@Smax, 1,
                        prob = model@localComm@Abundance[1:model@localComm@Smax])
            model@localComm <- birth(localComm)
        }
    }

    # # update equilib in params (if neccesary)
    # # stub
    #
    # # trim local community to max possible species
    # # JJ <- JJ[1:JiMax]

    # STUB: add random trait and pi info
    #comm$local_comm$Trait <- rnorm(length(comm$local_comm$JiMax))
    #comm$local_comm$pi = rlnorm(length(comm$local_comm$JiMax))

    # # add nstep to params
    # params$nstep <- nstep
    #
    # out <- list(local_comm = localComm, meta_comm = JJm, params = params)
    # class(out) <- 'roleComm'

    return(comm)
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

# testing
# params <- list(species_meta = 100,
#                individuals_meta = 10000,
#                individuals_local = 1000,
#                dispersal_prob = 0.1,
#                speciation_local = 0.01)
# testSim <- roleSimPlay(params, nstep = 100, nout = 5)
# testSim
c(-1:-30)
