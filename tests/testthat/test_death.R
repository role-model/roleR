library(testthat)
library(roleR)
context('death functions work')

# works
test_that('deathComm works', {
    l <- new(localCommCpp, rep(1,10),rep(1,10),10,rep(1,10))
    l$death(0)
    expect_equal(l$abundance_indv[1], 0)
})

# works
test_that('deathPhylo works', {
    p <- TreeSim::sim.bd.taxa(15, numbsim = 1,
                                lambda = 1,
                                mu = 0.8)[[1]]
    p <- apeToPhyloCpp(p)
    
    # 0 call here aligns with 1 index in alive 
    p$death(0)

    expect_false(p$alive[1])
})

#works - but double check comp filtering as we move forward
test_that('deathRole works', {
    r <- initSim()
    #r$local <- new(localCommCpp, rep(1,10),rep(1,10),10,rep(1,10))

    #prevtips <- sum(r$phylo$alive)

    r$death() 
    
    #verify that death occurred (there is a 0 abundance) 
    expect_true(is.element(0,r$local$abundance_indv))
})
