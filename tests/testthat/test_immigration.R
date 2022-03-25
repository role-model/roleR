library(testthat)
context('immigration functions work')

# works
test_that('immigrationLocal works', {
  
    # create a dummy localComm
    l <- new(localCommCpp, rep(1,10),rep(1,10),10,rep(1,10))
    
    # get a meta object to feed into imm function
    m <- scratchModel()
    meta <- m$meta
    
    # set an indv to dead
    l$abundance_indv[2] <- 0
    
    # call immigration on species 0, the immigrating individual placed at index 1
    l$immigration(0,1,m) 
    
    # check if abundance indv is back to 1
    expect_equal(l$abundance_indv[2], 1)
    
    # check if id of new species matches immigration species
    expect_equal(l$species_ids[2], 0)
    
    # check if local species abundance has been increased by 1
    expect_equal(l$abundance_sp[1], 2)
})

# works but uncertain
test_that('immigrationRole works', {
    r <- scratchModel()
    r$local <- new(localCommCpp, rep(1,10),rep(1,10),10,rep(1,10))
    
    # simulate death of indv at index 9 
    r$local$abundance_indv[9] = 0
    
    # call imm on death_index
    r$immigration(9)
  
    # there should now be an individual at that index
    expect_equal(r$local$abundance_indv[10],1)
    
    # the phylo tip should be alive again (I think???)
    expect_true(r$phylo$alive[10] == TRUE)
    
    # run it many times and see if it crashes
    r <- scratchModel()
    for(i in 1:100){
      d <- sample.int(15,1)
      r$immigration(d)
    }
    
})

sourceCpp("src/roleModelCpp.cpp")
sourceCpp("src/modules.cpp")
