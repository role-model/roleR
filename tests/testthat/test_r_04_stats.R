test_that("getSumStats defaults return several metrics", {
    
    p <- roleParams(niter = 200, speciation_local = .5) 
    # create a test model and get data from a step and the params
    m <- roleModel(p)
    m <- runRole(m)   
    stats <- getSumStats(m)
    expect_true(ncol(stats) > 1)
})
