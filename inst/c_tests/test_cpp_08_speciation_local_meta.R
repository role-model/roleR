
test_that("update_speciation_local_meta sets the new trait", {
})

test_that("update_speciation_local_meta sets the indv that was birthed or immigrated to the new species", {
    # create a test model and get data from a step and the params
    m <- quickModel()
    d <- m@modelSteps[[1]]
    p <- m@params
    
    # duplicate the data for comparison after the method
    d_start <- rlang::duplicate(d)

    #  update_speciation_local_meta, making the individual at dead_index 0 the new species
    d <- cFun(type="data",fun_name="update_speciation_local_meta",data=d,params=p,i=1,
                dead_index=0, dispersed_this_iter = TRUE,speciation_sp=0)

    # compare these two - the species of the new individual should equal the number of tips (as it is the newest tip)
    expect_true(d@localComm@indSpecies[1] == d@phylo@n)
})

test_that("update_speciation_local_meta updates the species last origin step", {
    # create a test model and get data from a step and the params
    m <- quickModel()
    d <- m@modelSteps[[1]]
    p <- m@params
    
    # duplicate the data for comparison after the method
    d_start <- rlang::duplicate(d)
    
    #  update_speciation_local_meta, making the individual at dead_index 0 the new species
    d <- cFun(type="data",fun_name="update_speciation_local_meta",data=d,params=p,i=8,
                dead_index=0, dispersed_this_iter = TRUE,speciation_sp=0)
    
    # since there were 10 tips before, the spLastOriginStep of the new
    #   species 11 should be 8, the iter we specified
    expect_true(d@localComm@spLastOriginStep[d_start@phylo@n + 1] == 8)
})
