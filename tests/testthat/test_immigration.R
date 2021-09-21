library(testthat)
context('immigration functions work')

test_that('immigrationLocal works', {
    c <- localComm(1:10,matrix(1:100,nrow = 10, ncol = 10),1:10,10)
    c <- immigration(c,1)
    expect_equal(c@abundance[1], 2)
})

test_that('immigrationRole works', {
    r <- .initSim(NULL)
    r@localComm <- localComm(1:10,matrix(1:100,nrow = 10, ncol = 10),1:10,10)
    r <- immigration(r)
    expect_equal(r@localComm@abundance[1],2)
})

