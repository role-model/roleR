
# WIP
test_that("update trait diffs squared changes the row of the dead indv", {
    
    # create a test model and get data from a step and the params
    m <- quickModel()
    d <- m@modelSteps[[2]]
    p <- m@params
    
    # add a new trait value at index 0 
    d@localComm@indTrait[1] <- 0.5

    # update trait diffs squared
    out <- cFun(type="data",fun_name="update_trait_diffs_sq",data=d,params=p,i=1, dead_index=0, parent_indv=0)
    # currently data does not contain trait diffs squared, so accessing this is difficult - TBD
    
    # compare these two - the trait diffs at the new trait value should be changed
})

test_that("update env filter probs does nothing if neutral", {
})

test_that("update env filter probs recalculates whole vector if env_sigma has changed this iter", {
})

test_that("update env filter probs sets one value of vector if env_sigma has NOT changed", {
})
