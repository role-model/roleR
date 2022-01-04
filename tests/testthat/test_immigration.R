library(testthat)
context('immigration functions work')

# works
test_that('immigrationLocal works', {
    l <- new(localCommCpp, rep(1,10),rep(1,10),10,rep(1,10))
    s <- initSim()
    m <- s$meta
    l$immigration(0,m) 
    expect_equal(l$abundance_indv[11], 1)
    expect_equal(l$species_ids[11], 0)
})

# works
test_that('immigrationRole works', {
    r <- initSim()
    r$local <- new(localCommCpp, rep(1,10),rep(1,10),10,rep(1,10))
    
    # simulate death of indv at index 9 
    r$local$abundance_indv[9] = 0
    # call imm on death_index
    r$immigration(9)
  
    # there should now be an individual at that index
    expect_equal(r$local$abundance_indv[10],1)
})

sourceCpp("src/roleModelCpp.cpp")
sourceCpp("src/modules.cpp")
