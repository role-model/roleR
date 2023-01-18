
test_that("birth sets the new indv to the parent species", {
    
    # create a test model and get data from a step and the params
    m <- quickModel()
    d <- m@modelSteps[[2]]
    p <- m@params
    
    # set the species of individual at index 1 to 1000
    d@localComm@indSpecies[1] <- 1000
    
    # duplicate the data for comparison after the call_birth method
    d_start <- rlang::duplicate(d)
    
    # call_birth using the params from iter = 1 on the dead_index at 0, with the parent at index 3
    # indices start from 0 
    out <- cFun(type="data",fun_name="call_birth",data=d,params=p,i=1, dead_index=1, parent_indv=0)
    
    # compare these two - the species of the new individual should be the same as the parent
    expect_true((out@localComm@indSpecies + 1)[2] == 1000)
})

test_that("birth adds 1 to the correct spot in the species matrix", {
    # create a test model and get data from a step and the params
    m <- quickModel()
    d <- m@modelSteps[[2]]
    p <- m@params
    
    # duplicate the data for comparison after the call_birth method
    d_start <- rlang::duplicate(d)
    
    # call_birth using the params from iter = 1 on the dead_index at 1, with the parent at index 0
    # indices start from 0 
    out <- cFun(type="data",fun_name="call_birth",data=d,params=p,i=1, dead_index=1, parent_indv=0)
    
    # get the abundance of the birthed species before the birth
    init_abundance <- d_start@localComm@spAbund[out@localComm@indSpecies[1] + 1] 
    # get after the birth
    out_abundance <- out@localComm@spAbund[out@localComm@indSpecies[1] + 1]

    # compare these two - the final abundance should be the start abundance plus 1 
    expect_true(out_abundance == init_abundance + 1)
})

test_that("birth adds a new trait value different from the last", {
    
    # create a test model and get data from a step and the params
    m <- quickModel()
    d <- m@modelSteps[[2]]
    p <- m@params
    
    # duplicate the data for comparison after the call_birth method
    d_start <- rlang::duplicate(d)
    
    # call_birth using the params from iter = 1 on the dead_index at 1, with the parent at index 0
    # indices start from 0 
    out <- cFun(type="data",fun_name="call_birth",data=d,params=p,i=1, dead_index=1, parent_indv=0)
    
    start_trait <- d_start@localComm@indTrait[2]
    end_trait <- out@localComm@indTrait[2]
    
    # compare these two
    expect_false(start_trait == end_trait)
})

# make sure it gives you a number in the allowed range 
test_that("birth calculates the new trait as expected", {
    # ie make sure its 0-1
})