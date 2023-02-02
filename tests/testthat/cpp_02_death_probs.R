
test_that("env filtering outputs a vector of the right length", {
    
    # create a test model and get data from a step and the params
    m <- quickModel()
    d <- m@modelSteps[[2]]
    p <- m@params
    
    # duplicate the data for comparison after the call_birth method
    d_start <- rlang::duplicate(d)
    
    # call_birth using the params from iter = 1 on the dead_index at 0, with the parent at index 3
    # indices start from 0 
    out <- cFun(type="vector",fun_name="get_filtering_death_probs",data=d,params=p,i=1)
    expect_true(m@params@niter == length(out))

    # try copying one liner into R
    # check edge cases 
    # make sure if you give it a negative number it doesnt work 
    # make sure vector is the right shape
    # wont be an opportunity to mess things up as we will check things in R
    
    # what is the expected behavior of this vector? can test using different params
})

test_that("env filtering probs is all the same if env_sigma param is 0 (no filtering)", {
    
    # create a test model and get data from a step and the params
    m <- quickModel()
    d <- m@modelSteps[[2]]
    p <- m@params
    
    # set param to 0
    p@env_sigma <- 0
    
    # get filtering death probs 
    out <- cFun(type="vector",fun_name="get_filtering_death_probs",data=d,params=p,i=0)
    
    expect_true(length(unique(out)) == 1)
})

test_that("env filtering probs is NOT all the same if env_sigma param is nonzero (filtering)", {
    
    # create a test model and get data from a step and the params
    m <- quickModel()
    d <- m@modelSteps[[2]]
    p <- m@params
    
    # set param to 0
    p@env_sigma <- 1
    
    # get filtering death probs 
    out <- cFun(type="vector",fun_name="get_filtering_death_probs",data=d,params=p,i=0)

    expect_false(length(unique(out)) == 1)
})