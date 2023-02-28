
trimPhyEdges <- function(e){
    return(e[e[,1] > -1,])
}

#WIP
test_that("update_speciation_phylo adds a new internal node", {
    # create a test model and get data from a step and the params
    m <- quickModel()
    p <- m@params
    d <- m@modelSteps[[1]]
    # duplicate the data for comparison after the call_birth method
    d_start <- rlang::duplicate(d)
    
    #  update_speciation_local_meta, making the individual at dead_index 0 the new species 
    out <- cFun(type="data",fun_name="update_speciation_phylo",data=d,params=p,i=1, 
                dead_index=0, dispersed_this_iter = TRUE,speciation_sp=0)
    
    trimPhyEdges(d_start@phylo@e)
    trimPhyEdges(out@phylo@e)
})

test_that("update_speciation_phylo adds new tips", {
})

test_that("update_speciation_phylo updates ancestry of internal nodes", {
})

test_that("update_speciation_phylo sets new edge length and augments old edge lengths", {
})

test_that("update_speciation_phylo DOESN'T update internal edges", {
})

test_that("update_speciation_phylo DOESN'T update extinct tips", {
})

test_that("update_speciation_phylo results in a phylo where the branches of all alive species lead to the present", {
})

test_that("update_speciation_phylo updates the alive vector", {
})

test_that("update_speciation_phylo increments the number of tips", {
})

#WIP
test_that("update_speciation_phylo creates a rolePhylo object that can be coerced to an ape object", {
    # create a test model and get data from a step and the params
    m <- quickModel()
    p <- m@params
    d <- m@modelSteps[[1]]
    # duplicate the data for comparison after the call_birth method
    d_start <- rlang::duplicate(d)
    
    # initial state tree works
    expect_error(plot(as(d_start@phylo, 'phylo')), NA)
    
    #  update_speciation_local_meta, making the individual at dead_index 0 the new species 
    d <- cFun(type="data",fun_name="update_speciation_phylo",data=d, params=p,i=1, 
                dead_index=0, dispersed_this_iter = TRUE,speciation_sp=0)
    
    trimPhyEdges(cbind(d_start@phylo@e,d@phylo@e))
    cbind(trimPhyEdges(d@phylo@e),trimPhyEdges(d_start@phylo@e))
    
    # tree works after one speciation
    expect_error(plot(as(d@phylo, 'phylo')), NA)
})

#WIP
test_that("update_speciation_phylo creates the expected phylogeny given a trivial 3 species tree and deterministic inputs", {
    p <- roleParams(individuals_local = 100, individuals_meta = 1000,
                    species_meta = 3, speciation_local = 1, 
                    speciation_meta = 0.5, extinction_meta = 0.05, env_sigma = 0.5,
                    trait_sigma=1,comp_sigma = 0.5, dispersal_prob = 0.1, mutation_rate = 0.01,
                    equilib_escape = 1, num_basepairs = 250,
                    init_type = 'oceanic_island', niter = 5, niterTimestep = 2)
    m <-  runRoLE(roleModel(p))
    p <- m@params
    d <- m@modelSteps[[1]]
    d1 <- m@modelSteps[[1]]
    d2 <- m@modelSteps[[2]] 
    d3 <- m@modelSteps[[3]]
    
    # duplicate the data for comparison after the call_birth method
    d_start <- rlang::duplicate(d)
    
    #  update_speciation_local_meta, making the individual at dead_index 0 the new species 
    # tree aint conformed unless call on sp 3
    d <- cFun(type="data",fun_name="update_speciation_phylo",data=d,params=p,i=1, 
                dead_index=0,speciation_sp=3)
    
    # plot start tree
    plot(as(d_start@phylo, 'phylo'))
    plot(as(d3@phylo,'phylo'))
    plot(as(d@phylo,'phylo'))
    
    # edges
    trimPhyEdges(d_start@phylo@e)
    trimPhyEdges(d@phylo@e)
    trimPhyEdges(d2@phylo@e)
    trimPhyEdges(d3@phylo@e)
    
    # compare n tips
    d_start@phylo@n
    d@phylo@n
    d2@phylo@n
    d3@phylo@n
    
    # tipnames
    d_start@phylo@tipNames[1:5]
    d@phylo@tipNames[1:5]
    d2@phylo@tipNames[1:5]
    
    # lengths
    d_start@phylo@l[1:10]
    d@phylo@l[1:10]
    d2@phylo@tipNames[1:5]
})
    