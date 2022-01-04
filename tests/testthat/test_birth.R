library(testthat)
library(roleR)

context('birth functions work')

# works
test_that('birthLocal works', {
    l <- new(localCommCpp, rep(1,10),rep(1,10),10,rep(1,10))
    l$birth(0)
    expect_equal(l$abundance_indv[11], 1)
})

# works
test_that('birthRole works', {
    r <- initSim()
    r$local <- new(localCommCpp, rep(1,10),rep(1,10),10,rep(1,10))
    r$birth()
    expect_equal(r$local$abundance_indv[11], 1)
})
