# simulate stochastic lotka-voltera

lvsim <- function(x0, params, nreps) {
    xt <-  matrix(NA, nrow = nreps, ncol = 2)
    xt[1, ] <- x0

    for(b in 2:nreps) {
        if(all(x0 <= 0)) break
        print(b)

        # rates
        births <- c(params$la1,
                    params$la2) * x0
        deaths <- c(params$mu1 * x0[1] + params$a12 * x0[2],
                    params$mu2 * x0[2] + params$a21 * x0[1]) * x0

        # which event happened
        e <- sample(1:4, 1, prob = c(births, deaths))

        # update
        if(e <= 2) {
            # brith
            i <- e
            x0[i] <- x0[i] + 1
        } else {
            # death
            i <- e - 2
            x0[i] <- x0[i] - 1
        }

        xt[b, ] <- x0
    }

    return(xt)
}

p <- list(la1 = 1,
          mu1 = 0.05,
          la2 = 0.8,
          mu2 = 0.05,
          a12 = 0.01,
          a21 = 0.01)

nn <- lvsim(x0 = c(10, 10), params = p, nreps = 1000)

matplot(nn, type = 'l', lty = 1, lwd = 2)

imax <- sum(!is.na(nn[, 1]))
colz <- viridis::viridis(imax)

plot(nn, type = 'n')

for(i in 2:imax) {
    segments(x0 = nn[i - 1, 1], x1 = nn[i, 1],
             y0 = nn[i - 1, 2], y1 = nn[i, 2],
             col = colz[i])
}







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



