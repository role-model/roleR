test_that("rcpp competition matrix calculation works", {
    # trait data, rows are individuals
    trt <- cbind(1:3, 3:5, 5:7)
    
    # competitive distance matrix
    s <- 1.5
    d <- roleR:::compMatCalcTest(x = trt, sigC = s)
    
    # distances calculated by hand
    dh <- dist(trt) |> as.matrix()
    dh <- exp(-(dh / s)^2)[2:3, 1]
    names(dh) <- NULL
    
    expect_equivalent(list("diag_equal" = diag(d), 
                           "up_low_equal" = d[lower.tri(d)], 
                           "dist_correct" = d[2:3, 1]), 
                      list("diag_equal" = rep(0, 3), 
                           "up_low_equal" = d[upper.tri(d)], 
                           "dist_correct" = dh))
})


test_that("rcpp competition matrix calculation works when all traits equal", {
    # trait data, rows are individuals
    trt <- matrix(-3, nrow = 3, ncol = 1)
    
    # competitive distance matrix
    s <- 1.5
    d <- roleR:::compMatCalcTest(x = trt, sigC = s)
    
    # distances calculated by hand
    dh <- dist(trt) |> as.matrix()
    dh <- exp(-(dh / s)^2)[2:3, 1]
    names(dh) <- NULL
    
    expect_equivalent(list("diag_equal" = diag(d), 
                           "up_low_equal" = d[lower.tri(d)], 
                           "dist_correct" = d[2:3, 1]), 
                      list("diag_equal" = rep(0, 3), 
                           "up_low_equal" = d[upper.tri(d)], 
                           "dist_correct" = dh))
})

