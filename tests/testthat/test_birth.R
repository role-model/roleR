library(testthat)
context('birth functions work')

test_that('birthLocal works', {
    l <- new(localCommCpp, rep(1:10), matrix(1:100, nrow = 10), 10, rep(1:10))
    l$birth(0)
    expect_equal(l$abundance[1], 2)
})

test_that('birthRole works', {
    r <- initSim()
    r$local <- new(localCommCpp, rep(1:10), matrix(1:100, nrow = 10), 10, rep(1:10))
    r$birth()
    
    expect_true(is.element(2,r$local$abundance))
})
