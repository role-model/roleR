

test_that("sample_zero_to_x draws randomly with the expected median", {
    
    # draw 100 times from indices 0-4
    out <- c()
    for(i in 1:100){
        sampled <- cFun(type="int",fun_name="sample_zero_to_x",probs=c(0),x=5) 
        out <- c(out,sampled)
    }

    # median should be 2 
    expect_true(2 == median(out))
})

test_that("sample_index_using_probs draws properly when probs are 0 or 1", {
    
    # guaranteed probs where the draw should always be index 2 
    probs <- c(0,0,1)
    sampled <- cFun(type="int",fun_name="sample_index_using_probs",probs=probs,x=0) 
    
    # draw should be 2 
    expect_true(2 == sampled)
})