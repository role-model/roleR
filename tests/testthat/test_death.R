library(testthat)
context('death functions work')

test_that('deathComm works', {
    l <- new(localCommCpp, rep(1:10), matrix(1:100, nrow = 10), 10, rep(1:10))
    l$death(0)
    expect_equal(l$abundance[1], 0)
})

test_that('deathPhylo works', {
    p <- TreeSim::sim.bd.taxa(15, numbsim = 1,
                                lambda = 1,
                                mu = 0.8)[[1]]
    p <- .apeToPhyloCpp(p)
    
    p$death(1)
    #unfinished
    expect_false(p$alive[1])
})

test_that('deathRole works', {
    r <- initSim()
    #set starting species (starts at 100) to have abundance of 1
    index <- match(100,r$local$abundance)
    r$local$abundance[index] <- 1
    r$death()
    
    #verify that death occured (no 1s)
    expect_true(!is.element(1,r$local$abundance))
    #unfinished
    #verify that phylodeath occurred, i.e. tip has been set to dead
    #below does not work currently due to misaligned indices, todo figure this out
    expect_false(r$phylo$alive[index])
})
