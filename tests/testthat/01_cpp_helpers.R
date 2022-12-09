test_that("sample_zero_to_x draws randomly with the expected median", {
    
    out <- c()
    for(i in 1:100){
        sampled <- cFun(type="int",fun_name="sample_zero_to_x",x=5) 
        out <- c(out,sampled)
    }
    
    expect_equal(3, median(out))
})

test_that("sample_index_using_probs draws properly when probs are 0 or 1", {
    
    probs <- c(0,0,1)
    sampled <- cFun(type="int",fun_name="sample_index_using_probs",probs=probs) 
    
    expect_equal(2, sampled)
})
