library(testthat)
library(roleR)

context('birth functions work')

# works
test_that('birthLocal works', {
  
    # create dummy local comm
    l <- new(localCommCpp, rep(1,10),rep(1,10),10,rep(1,10))
    l$traits[2] <- 2
    
    # set an individual to dead
    l$abundance_indv[1] <- 0
    
    # look at some values 
    l$abundance_indv
    l$abundance_sp
    l$species_ids
    l$traits
    
    # call birth on individual 1, the new indv replacing the indv at index 0
    l$birth(1,0)
    
    # look at the new values
    l$abundance_indv
    l$abundance_sp
    l$species_ids
    l$traits
    
    # the individual at the dead index should now be alive
    expect_equal(l$abundance_indv[1], 1)
    
    # the species id at the dead index should now be the species of the birth
    expect_equal(l$species_ids[1], 1)
    
    # the trait at the dead index should now be different
    expect_equal(l$traits[1], 2)
    
    # the species of the individual that gave birth (1) should be incremented by 1
    expect_equal(l$abundance_sp[2], 2)
})

# works
test_that('birthRole works', {
  
    # create a dummy params and make a sim
    params <- roleParams(nrun=1,niter=10,niterTimestep=10,defaults=TRUE)
    cparams <- stretchAndSampleParams(params)
    parlist <- cparams@values[[1]]
    sim <- initSim(parlist,type="bridge_island")
    
    # assign a dummy local comm
    sim$local <- new(localCommCpp, rep(1,10),rep(1,10),10,rep(1,10))
    
    # call birth on dead index 1
    # pretty much just does local birth
    sim$birth(1)
    expect_equal(sim$local$abundance_indv[1], 1)
    
    # run it many times and see if it crashes
    r <- scratchModel()
    for(i in 1:100){
      d <- sample.int(119,1) # verify that model$birth(i) maxes at i = J-1 
      r$birth(d)
    }
})
