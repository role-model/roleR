library(testthat)

context('speciation functions work')

# works
test_that('specLocal works', {
  
    r <- scratchModel()
    
    l <- new(localCommCpp, rep(1,10),rep(1,10),10,rep(1,10))
    
    #l$species_ids
    #l$abundance_indv
    #l$abundance_sp
    #l$Smax
    
    # call speciation of species i replacing indv at index 9  
    l$speciation(1,9,r$params$trait_sigma[1])
    
    #l$species_ids
    #l$abundance_indv
    #l$abundance_sp
    
    # species of new indv should be i 
    expect_equal(l$species_ids[10], 10)
    
    # Smax should be 11 
    expect_equal(l$Smax,11)
})

# works
test_that('specPhylo works', {
    pp <- TreeSim::sim.bd.taxa(15, numbsim = 1,
                              lambda = 1,
                              mu = 0.8, complete = FALSE)[[1]]
    p <- apeToPhyloCpp(pp)
    #p$alive
    #p$n
    length(p$alive)
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
    
    # run it many times and see if it crashes
    r <- scratchModel()
    for(i in 1:100){
      d <- sample.int(120,1) - 1 # verify that model$birth(i) maxes at i = J-1 
      r$speciation(d)
    }
})

