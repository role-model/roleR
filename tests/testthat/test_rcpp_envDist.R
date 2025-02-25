test_that("rcpp enviro distance matrix calculation works", {
    # trait data, rows are individuals
    trt <- cbind(1:3, 3:5, 5:7) / 10
    
    # competitive distance matrix
    s <- 2
    eOpt <- matrix(c(0, 0, 1), nrow = 1)
    d <- roleR:::envDistCalcTest(x = trt, envOptim = eOpt, sigE = s)
    
    # distances calculated by hand
    dh <- dist(rbind(eOpt, trt)) |> as.matrix()
    dh <- exp(-(dh / s)^2)[2:4, 1, drop = FALSE]
    rownames(dh) <- colnames(dh) <- NULL
    
    expect_equivalent(d, dh)
})
