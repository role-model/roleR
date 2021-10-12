library(testthat)
context('speciation functions work')

test_that('specLocal works', {
    r <- .initSim(NULL)
    c <- localComm(rep(1,10),matrix(1:100,nrow = 10, ncol = 10),1:10,10)
    c <- speciation(c, 1, r@params)
    expect_equal(c@Smax, 11)
    expect_equal(c@abundance[11],1)
})

test_that('specPhylo works', {
    p <- TreeSim::sim.bd.taxa(15, numbsim = 1,
                              lambda = 1,
                              mu = 0.8, complete = FALSE)[[1]]
    p <- as(p, "rolePhylo")
    p <- speciation(p,1)
    expect_true(p@alive[16])
})

test_that('specRole works', {
    r <- .initSim(NULL)

    prevSumAlive <- sum(r@phylo@alive)
    prevSmax <- r@metaComm@Smax

    r <- speciation(r)
    expect_true(sum(r@phylo@alive) == prevSumAlive + 1)
    expect_true(sum(r@metaComm@Smax) == prevSmax + 1)


})

