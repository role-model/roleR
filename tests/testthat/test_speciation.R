library(testthat)
context('speciation functions work')

test_that('specComm works', {
    c <- localComm(rep(1,10),matrix(1:100,nrow = 10, ncol = 10),1:10,10)
    c <- speciation(c)
    expect_equal(c@Smax, 11)
    expect_equal(c@abundance[11],1)
})

test_that('specLocal works', {
    c <- localComm(rep(1,10),matrix(1:100,nrow = 10, ncol = 10),1:10,10)
    c <- speciation(c)
    expect_equal(c@Smax, 11)
    expect_equal(c@abundance[11],1)
    #unfinished
})

test_that('specPhylo works', {
    p <- TreeSim::sim.bd.taxa(15, numbsim = 1,
                              lambda = 1,
                              mu = 0.8)[[1]]
    p <- as(p, "rolePhylo")
    p <- speciation(p,1)
    #unfinished
    #expect_true(p@alive[16])
})


#unfinished
test_that('specRole works', {
    r <- .initSim(NULL)
    r@localComm <- localComm(rep(1,10),matrix(1:100,nrow = 10, ncol = 10),1:10,10)
    source("R/speciation.R")
    r <- speciation(r)
    #unfinished
})

