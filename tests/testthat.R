library(testthat)
library(roleR)

# function used by multiple tests
# take a run role model object and return a time-series-like
# matrix where rows are time steps and there is a column for each spp
# cells are proportional abundances

ats <- function(r) {
    nspp <- r@modelSteps[[length(r@modelSteps)]]@phylo@n
    
    res <- lapply(r@modelSteps, function(x) {
        sapply(1:nspp, function(i) {
            mean(x@localComm@indSpecies == i)
        })
    })
    
    return(do.call(rbind, res))
}

set.seed(1)
test_check('roleR')
set.seed(NULL)
