#' simulate IBM lotka-voltera
#' 
#' @param x0 initial condition, vector of spp IDs (1s and 2s)
#' @param t0 initial trait value of each individual
#' @param params a named list of the model params:
#'     * `de` 2-long vector of death rates due to environment
#'     * `dc` 2-long vector of death rates due to competition
#'     * `se` sigma for environment death rate kernel 
#'     * `sc` sigma for competition death rate kernel 
#'     * `that` "t hat" the optimal trait value
#'     * `sbm` rate of brownian trait variation
#' @param nreps single int, the number of timesteps
#' @return matrix of nrow equal to `nreps` and ncol equal to 2 (one for each spp)


lvibm <- function(x0, t0, params, nreps) {
    # matrix to store abundances
    xt <-  matrix(NA, nrow = nreps, ncol = 2)
    xt[1, ] <- c(sum(x0 == 1), sum(x0 == 2))

    # vector of trait differences to optimum
    we <- (t0 - params$that)^2

    # matrix of sp-sp trait diffs
    wc <- as.matrix(dist(t0))^2

    for(b in 2:nreps) {
        if(all(x0 <= 0)) break
        print(b)

        # alive or dead
        alive <- x0 > 0
        ialive <- which(alive)
        xalive <- x0[alive]
        talive <- t0[alive]
        wealive <- we[alive]
        wcalive <- wc[alive, alive, drop = FALSE]

        # birth rate
        births <- 1

        # relative death rates
        # browser()
        eComp <- params$de[xalive] * (1 - exp(-params$se * wealive))
        # wtf <- try(rowSums(exp(-wcalive / params$sc)))
        # if(class(wtf) == 'try-error') browser()
        cComp <- params$dc[xalive] *
            as.vector(rowSums(exp(-wcalive / params$sc)))
        deaths <- eComp + cComp

        # birth or death
        e <- sample(1:2, 1, prob = c(births, sum(deaths)))

        # update
        if(e == 1) {
            # birth
            # who gives birth
            j <- sample(length(ialive), 1)
            iparent <- ialive[j]
            inew <- which(!alive)[1]
            x0[inew] <- x0[iparent]
            t0[inew] <- t0[iparent] + rnorm(1, 0, sd = sqrt(params$sbm))

            # update trait distance objects
            we[inew] <- (t0[inew] - p$that)^2
            wc[inew, ] <- (t0[inew] - t0)^2
            wc[, inew] <- (t0[inew] - t0)^2
        } else {
            # death
            # who died
            j <- sample(length(ialive), 1, prob = deaths)
            idie <- ialive[j]
            x0[idie] <- 0
        }

        # update output
        # browser()
        xt[b, ] <- c(sum(x0 == 1), sum(x0 == 2))
    }

    return(xt)
}



p <- list(de = c(0, 0),
          dc = c(0.02, 0.02),
          se = 2,
          sc = 0.1,
          that = 0,
          sbm = 0.1)

x <- c(rep(1:2, each = 10), rep(0, 980))
trt <- c(rnorm(10, -2), rnorm(10, 2), rep(0, 980))

mm <- lvibm(x0 = x, t0 = trt, params = p, nreps = 5000)

matplot(mm, type = 'l', lty = 1)



