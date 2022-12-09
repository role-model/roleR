
test_that("birth sets the new indv to the parent species", {
    # birth
    d <- model@modelSteps[[2]]
    d@localComm@indSpecies[1] <- 1
    p <- model@params
    
    d1 <- rlang::duplicate(d)
    d@localComm@indSpecies
    
    out <- cFun(type="data",fun_name="call_birth",data=d,params=p,i=1, dead_index=0, parent_indv=3)
    
    out@localComm@indSpecies
    d1@localComm@indSpecies
})



test_that("birth sets the new indv to the parent species", {
    
})