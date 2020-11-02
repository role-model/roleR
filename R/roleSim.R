#' @title RoLE model simulation
#'
#' @description Simulate communities under the RoLE model
#'
#'
#' @param params list of parameters
#' @param nstep number of simulation steps to run
#' @param nsim number of simulations to run
#'
#' @details Stub
#'
#' @return An object of class \code{roleComm}
#' @export

roleSim <- function(params, nstep = NULL, nsim = 1) {
    # number of species in the metacommunity
    Sm <- params$species_meta

    # metacommunity SAD
    JJm <- pika::rfish(Sm, 0.01)
    JJm <- JJm / sum(JJm)

    # number of local individuals
    J <- params$individuals_local

    # vector of local species abundances
    JJ <- rep(0, Sm * 100)

    # initialize local species abundances with one species having all individuals
    JJ[sample(Sm, 1, prob = JJm)] <- J

    # counter keeping track of max number of possible species in local community
    JiMax <- Sm

    # further parameters
    m <- params$dispersal
    nu <- params$speciation_local

    nT <- nstep

    for(i in 1:nT) {
        # death
        dead <- sample(JiMax, 1, prob = JJ[1:JiMax])
        JJ[dead] <- JJ[dead] - 1


        if(runif(1) <= nu) {
            # speciation
            JJ[JiMax + 1] <- 1
            JiMax <- JiMax + 1
        } else if(runif(1) <= m) {
            # immigration
            imm <- sample(Sm, 1, prob = JJm)
            JJ[imm] <- JJ[imm] + 1
        } else {
            # local birth
            birth <- sample(JiMax, 1, prob = JJ[1:JiMax])
            JJ[birth] <- JJ[birth] + 1
        }
    }

    # update equilib in params (if neccesary)
    # stub

    # add nstep to params
    params$nstep <- nstep

    out <- list(local_comm = JJ, params = params)
    class(out) <- 'roleComm'

    return(out)

}
