library(testthat)
context('immigration functions work')

test_that('immigrationLocal works', {
    l <- new(localCommCpp, rep(1,10),rep(1,10),10,rep(1,10))
    s <- initSim()
    m <- s$meta
    l$immigration(0,m) 
    expect_equal(l$abundance_indv[11], 1)
    expect_equal(l$species_ids[11], 0)
})

test_that('immigrationRole works', {
    r <- initSim()
    r$local <- new(localCommCpp, rep(1:10), matrix(1:100, nrow = 10), 10, rep(1:10))
    r$immigration()
    
    expect_true(is.element(2,r$local$abundance))
    r <- .initSim(NULL)
    r@localComm <- localComm(1:10,matrix(1:100,nrow = 10, ncol = 10),1:10,10)
    r <- immigration(r)
    expect_equal(r@localComm@abundance[1],2)
})

sourceCpp("src/roleModelCpp.cpp")
sourceCpp("src/modules.cpp")
