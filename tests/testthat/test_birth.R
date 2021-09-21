library(testthat)
context('birth functions work')

test_that('birthLocal works', {
    c <- localComm(1:10,matrix(1:100,nrow = 10, ncol = 10),1:10,10)
    c <- birth(c,1)
    expect_equal(foo@abundance[1], 2)
})

#unfinished
test_that('birthRole works', {
    r <- .initSim(NULL)
    r@localComm <- localComm(rep(1,10),matrix(1:100,nrow = 10, ncol = 10),1:10,10)
    r <- birth(r)
    expect_true(is.element(2,r@localComm@abundance))
})

