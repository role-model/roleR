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
    init$params <- params

    # ----
    # loop over steps
    .iterSim(init, nstep)
}



#' @rdname roleSim
#' @export

roleSimPlay <- function(params, init = NULL, nstep = NULL, nout = NULL) {
    # ----
    # initialize simulation

    if(is.null(init)) {
        init <- .initSim(params)
    }

    # make sure `params` reflect the user-specified values (potentially
    # overwriting old values)
    init$params <- params

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
    # empty `roleComm` object
    out <- list(local_comm = list(Abundance = numeric(),
                                  Traits = numeric(),
                                  pi = numeric()),
                meta_comm = list(Abundance = numeric(),
                                 Traits = numeric()),
                phylo = list(),
                params = params)
    class(out) <- 'roleComm'

    # metacommunity SAD
    out$meta_comm$Abundance <- .lseriesFromSN(params$species_meta,
                                              params$individuals_meta)

    # initialize phylo
    out$phylo <- TreeSim::sim.bd.taxa(params$species_meta, numbsim = 1,
                                      lambda = params$speciation_meta,
                                      mu = params$extinction_meta)[[1]]

    # General not on book keeping: for vectors that we expect to grow (i.e.,
    # traits and local abundances) we make a vector 100 times longer than the
    # initial number of species in the metacommunity rather than needing to
    # augment those vectors

    # initialize traits
    out$Traits <- c(ape::rTraitCont(out$phylo, sigma = params$sigma_bm),
                    rep(NA, params$species_meta * 99))

    # vector of local species abundances
    out$local_comm$Abundance <- rep(0, params$species_meta * 100)

    # initialize local species abundances with one species having all
    # individuals
    out$local_comm$Abundance[sample(params$species_meta, 1,
                                    prob = out$meta_comm$Abundance)] <-
        params$individuals_local

    # counter keeping track of max number of possible species in local
    # community
    out$local_comm$JiMax <- params$species_meta

    return(out)
}


# ----
# function to iterate role simulation
#' @param comm an object of class `roleComm`
#' @param nstep number of steps to iterate

.iterSim <- function(comm, nstep) {
    # browser()
    for(i in 1:nstep) {
        # death
        dead <- sample(comm$local_comm$JiMax, 1,
                       prob = comm$local_comm$Abundance[1:comm$local_comm$JiMax])
        comm$local_comm$Abundance[dead] <- comm$local_comm$Abundance[dead] - 1

        if(runif(1) <= comm$params$speciation_local) {
            # speciation
            comm$local_comm$Abundance[comm$local_comm$JiMax + 1] <- 1
            comm$local_comm$JiMax <- comm$local_comm$JiMax + 1

            # need to update phylogeny
            # need to assign new trait
        } else if(runif(1) <= comm$params$dispersal_prob) {
            # immigration
            imm <- sample(comm$params$species_meta, 1,
                          prob = comm$meta_comm$Abundance)
            comm$local_comm$Abundance[imm] <- comm$local_comm$Abundance[imm] + 1
        } else {
            # local birth
            birth <-
                sample(comm$local_comm$JiMax, 1,
                       prob = comm$local_comm$Abundance[1:comm$local_comm$JiMax])
            comm$local_comm$Abundance[birth] <-
                comm$local_comm$Abundance[birth] + 1
        }
    }

    # # update equilib in params (if neccesary)
    # # stub
    #
    # # trim local community to max possible species
    # # JJ <- JJ[1:JiMax]

    # STUB: add random trait and pi info
    comm$local_comm$Trait <- rnorm(length(comm$local_comm$JiMax))
    comm$local_comm$pi = rlnorm(length(comm$local_comm$JiMax))

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
