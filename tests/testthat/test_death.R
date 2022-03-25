library(testthat)
library(roleR)
context('death functions work')

# works
test_that('deathComm works', {
  
    # create dummy local comm
    l <- new(localCommCpp, rep(1,10),rep(1,10),10,rep(1,10))
    # call death on indv 0
    l$death(0)
    # verify that individual is dead (note that 1 index in R equals 0 index in C++)
    expect_equal(l$abundance_indv[1], 0)
})

test_that('deathPhylo works', {
  
    # create test tree
    p <- TreeSim::sim.bd.taxa(15, numbsim = 1,
                                lambda = 1,
                                mu = 0.8, complete=FALSE)[[1]]
    p <- apeToRolePhylo(p)
    p <- rolePhyloToCpp(p)
    
    p <- rolePhyloFromCpp(p)
    
    p$alive
    
    # 0 call here aligns with 1 index in alive 
    # call death
    p$death(1)
    
    p$alive
    
    # verify dead
    expect_false(p$alive[1])
})

#works - but double check comp filtering as we move forward
test_that('deathRole works', {
  
    # create a dummy params and make a sim
    # replace with scratchModel()
    params <- roleParams(nrun=1,niter=10,niterTimestep=10,defaults=TRUE)
    cparams <- stretchAndSampleParams(params)
    parlist <- cparams@values[[1]]
    sim <- initSim(parlist,type="bridge_island")
    sim$print <- TRUE
    
    # assign a dummy local comm
    sim$local <- new(localCommCpp, rep(1,10),rep(1,10),10,rep(1,10))
  
    # look at values
    sim$local$abundance_indv
    sim$local$abundance_sp
    sum(sim$phylo$alive == FALSE)
    
    # death crashes it sometimes
    # call death, which samples for an index
    sim$death()
 
    sim$local$abundance_indv
    
    #verify that there is a 0 abundance in abundance_indv
    expect_true(is.element(0,sim$local$abundance_indv))

    #verify that there is a 0 abundance in abundance_sp (a species has been decremented)
    expect_true(is.element(0,sim$local$abundance_sp))
    
    #verify that there is a new FALSE in alive (a phylo death was triggered due to a species going extinct)
    expect_true(sum(sim$phylo$alive == FALSE) == 1)
    
    # run it many times and see if it crashes
    r <- scratchModel()
    for(i in 1:100){
      r$death()
    }
})
