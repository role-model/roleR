library(testthat)
context('speciation functions work')

test_that('specLocal works', {
    r <- initSim()
    l <- r$local
    l$speciation(1,r$params)
    expect_equal(l$Smax,16)
})

test_that('specPhylo works', {
    p <- TreeSim::sim.bd.taxa(15, numbsim = 1,
                              lambda = 1,
                              mu = 0.8, complete = FALSE)[[1]]
    p <- .apeToPhyloCpp(p)
    p$speciation(0)
    expect_true(p@alive[16])
})

test_that('specRole works', {
    r <- initSim()

    prevSumAlive <- sum(r$phylo$alive)
    prevSmax <- r$meta$Smax

    r$speciation()
    expect_true(sum(r@phylo@alive) == prevSumAlive + 1)
    expect_true(sum(r@metaComm@Smax) == prevSmax + 1)
})

