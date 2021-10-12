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
    #set starting species (starts at 100) to have abundance of 1
    index <- match(100,r@localComm@abundance)
    index
    r@localComm@abundance[index] <- 1
    r@localComm@abundance
    r@phylo@alive
    r <- death(r)
    r@phylo@alive
    #verify that death occured (no 1s)
    expect_true(!is.element(1,r@localComm@abundance))
    #verify that phylodeath occurred, i.e. tip has been set to dead
    #below does not work currently due to misaligned indices, todo figure this out
    expect_false(r@phylo@alive[index])
})
