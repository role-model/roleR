#' simulate stochastic lotka-voltera
#' 
#' @param x0 initial condition, should be 2-long int vector
#' @param params a named list of the model params:
#'     * `la1` birth rate for sp 1
#'     * `mu1` density dependent death rate for sp 1
#'     * `la2` birth rate for sp 2
#'     * `mu2` density dependent death rate for sp 2
#'     * `a12` density dependent effect of sp 2 on 1
#'     * `a21` density dependent effect of sp 1 on 2
#' @param nreps single int, the number of timesteps
#' @return matrix of nrow equal to `nreps` and ncol equal to 2 (one for each spp)


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


