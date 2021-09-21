library(testthat)
context('death functions work')

test_that('deathComm works', {
    c <- localComm(1:10,matrix(1:100,nrow = 10, ncol = 10),1:10,10)
    c <- death(c,1)
    expect_equal(c@abundance[1], 0)
})

test_that('deathPhylo works', {
    p <- TreeSim::sim.bd.taxa(15, numbsim = 1,
                                lambda = 1,
                                mu = 0.8)[[1]]
    p <- as(p, "rolePhylo")
    p <- death(p,1)
    expect_false(p@alive[1])
})

test_that('deathRole works', {
    r <- .initSim(NULL)
    r <- death(r)
    #unfinished
})

