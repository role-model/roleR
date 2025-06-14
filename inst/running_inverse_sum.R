# figure out running inverse sum ----

# function to update inverse sum
#' @param jdead index of dead *sp*
#' @param jborn index of born *sp*
#' @param s holds running inverse sum
#' @param lu indicates last time sp was updated
#' @param a holds running species abundance
#' @param step index of step we're on

updateIS <- function(jdead, jborn, s, lu, a, step) {
    if(jdead == jborn) {
        # death replaced by birth all in same sp, so no need update `a`
        # update from previous steps
        # if (last update step) == (step - 1), donʻt do anything
        # we achive this by `* 0`
        s[jborn] <- s[jborn] + (1 / a[jborn]) * (step - lu[jborn] - 1)
        
        # update inverse sum for the current step
        s[jborn] <- s[jborn] + (1 / a[jborn]) 
    } else {
        # update dead one from previous steps
        s[jdead] <- s[jdead] + (1 / a[jdead]) * (step - lu[jdead] - 1)
        
        # if born one had positive abund prior, update it from previous steps
        if(a[jborn] > 0) {
            s[jborn] <- s[jborn] + (1 / a[jborn]) * (step - lu[jborn] - 1)
        }
        
        # update abundances
        a[jborn] <- a[jborn] + 1
        a[jdead] <- a[jdead] - 1
        
        # update inv sum for the current step
        s[jborn] <- s[jborn] + (1 / a[jborn]) 
        
        if(a[jdead] > 0) { # only updated inv sum if not extirpated
            s[jdead] <- s[jdead] + (1 / a[jdead])
        } else { # if extirpated, set inv sum to 0
            s[jdead] <- 0
        }
    }
    
    return(list(s, a))
}

# simulation set-up
spAbund <- c(10, 0, 0, 0)
invSum <- 1 / spAbund
invSum[!is.finite(invSum)] <- 0
lastUpdate <- lastOrig <- c(0, 0, 0, 0)
nstep <- 11

spMat <- matrix(0, nrow = nstep + 1, ncol = length(spAbund))
rownames(spMat) <- 0:nstep
colnames(spMat) <- 1:4
spMat[1, ] <- spAbund

# store intermediate steps
kwriteOut <- c(1, 4, 8, 12)

outList <- vector(mode = "list", length = length(kwriteOut))
names(outList) <- paste("step", kwriteOut - 1, sep = "_")

outList[[1]] <- list(0, invSum, lastUpdate, lastOrig)

for(k in 2:(nstep + 1)) {
    # k is index for matrix
    # i is step index
    i <- k - 1
    
    # set up prob of death so only spp > 0 abund can experience death
    # but otherwise death prob is uniform (thus more likely to have
    # extirpation so we can test that scenario)
    dprob <- 1 * (spAbund > 0)
    
    # who dead, who born
    idead <- sample(4, 1, prob = dprob) 
    iborn <- sample(4, 1)
    
    if(spAbund[iborn] == 0) {
        lastOrig[iborn] <- i
    }
    
    u <- updateIS(idead, iborn, invSum, lastUpdate, spAbund, i)
    
    invSum <- u[[1]]
    spAbund <- u[[2]]
    
    lastUpdate[idead] <- i
    lastUpdate[iborn] <- i
    
    # if the we are writing out, make sure to update everyone not already 
    # updated on this step
    if(k %in% kwriteOut) {
        for(j in 1:4) { # j is a species index
            if(lastUpdate[j] != k - 1) {
                # take special care for 0 abundance becase when
                # `jdead == jborn`, `updateIS` method doesnʻt work for 0 abund
                if(spAbund[j] > 0) {
                    u <- updateIS(j, j, invSum, lastUpdate, spAbund, i)
                    
                    # only need to update `invSum` because `spAbund` 
                    # hasn't changed
                    invSum <- u[[1]] 
                } else {
                    invSum[j] <- 0
                }
                
                # update `lastUpdate`
                lastUpdate[j] <- i
                
                # write-out stuff
                outList[[paste("step", i, sep = "_")]] <- 
                    list(invSum = invSum, 
                         lastUpdate = lastUpdate, 
                         lastOrig = lastOrig)
            }
        }
    }
    
    spMat[k, ] <- spAbund
}


data.frame(writeout = (1:nrow(spMat) %in% kwriteOut), spMat)

for(i in 2:length(kwriteOut)) {
    pls <- apply(spMat[1:kwriteOut[i], ], 2, function(x) {
        cis <- 0
        for(j in x) {
            if(j == 0) {
                cis <- 0
            } else {
                cis <- cis + 1 / j
            }
        }
        
        cis
    })
    
    print(outList[[names(outList)[i]]]$invSum == pls)
}


# figure out if we have to keep track of harm mean of proportion ----

hmean <- function(x) {
    length(x) / sum(1 / x)
}

npop <- 50
mJ <- 100
px <- 0.5

J <- rpois(10, mJ)
x <- rpois(10, px * mJ)

all(x < J)

hmean(x / J) - hmean(x) / hmean(J)


# the harm mean of the proportion and hmean(x) / hmean(J) are close enough
