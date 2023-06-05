
test_that("dispersal sets the species as the parent species from meta", {
    
    # create a test model and get data from a step and the params
    m <- quickModel()
    d <- m@modelSteps[[2]]
    p <- m@params
    
    # should set the species at index 0 of indSpecies to 9, the meta species 
    out <- cFun(type="data",fun_name="call_dispersal",data=d,params=p,i=1, dead_index=0, parent_indv=9)
    
    # compare these two - the species of the new individual should be the same as the parent
    expect_true((out@localComm@indSpecies)[1] == 9)
})


test_that("dispersal increments abundance", {
    
    # create a test model and get data from a step and the params
    m <- quickModel()
    d <- m@modelSteps[[2]]
    p <- m@params
    
    # should set the species at index 0 of indSpecies to 9, the meta species 
    out <- cFun(type="data",fun_name="call_dispersal",data=d,params=p,i=1, dead_index=0, parent_indv=9)
    out@localComm@spAbund[10]
    
    # compare these two - the species of the new individual should be the same as the parent
    expect_true(out@localComm@spAbund[10] == 1)
})

test_that("dispersal updates the last time of the species origin", {
    
    # create a test model and get data from a step and the params
    m <- quickModel()
    d <- m@modelSteps[[2]]
    p <- m@params
    
    # call dispersal on species index 9 in the metacommunity at timestep 1
    # this should cause the last time of origin of species 9 to be 1 
    out <- cFun(type="data",fun_name="call_dispersal",data=d,params=p,i=1, dead_index=0, parent_indv=9)
    expect_true(out@localComm@spLastOriginStep[10] == 1)
})
