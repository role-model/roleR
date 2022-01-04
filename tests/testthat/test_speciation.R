library(testthat)

context('speciation functions work')

# works
test_that('specLocal works', {
    r <- initSim()
    l <- r$local
    
    # call speciation of species i replacing indv at index 9  
    l$speciation(1,9,r$params)
    
    # species of new indv should be i 
    l$species_ids[10]
    
    # Smax should be 16 
    expect_equal(l$Smax,16)
})

# works
test_that('specPhylo works', {
    pp <- TreeSim::sim.bd.taxa(15, numbsim = 1,
                              lambda = 1,
                              mu = 0.8, complete = FALSE)[[1]]
    p <- apeToPhyloCpp(pp)
    p$n
    p$speciation(1)
    expect_true(p$alive[15])
})

# works
test_that('specRole works', {
    r <- initSim()
    prevSumAlive <- sum(r$phylo$alive)
    prevSmax <- r$meta$Smax
    r$speciation(9)
    
    expect_true(sum(r$phylo$alive) == prevSumAlive + 1)
    expect_true(sum(r$local$Smax) == prevSmax + 1)
})

