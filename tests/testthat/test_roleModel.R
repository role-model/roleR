#test_roleModel.R
library(testthat)
library(roleR)
context('several model runs work')

# works
test_that('model of 100 iters runs without crash with filled timeseries', {
  model <- roleR::dummyModel(R=FALSE,run=TRUE,fill_ts=TRUE,niters=100)
})

# works
test_that('model of 100 iters runs without crash', {
  model <- roleR::dummyModel(R=FALSE,run=TRUE,fill_ts=FALSE,niters=100)
})

# works
test_that('model of 1000 iters runs without crash', {
  model <- dummyModel(R=FALSE,run=TRUE,niters=1000)
})
