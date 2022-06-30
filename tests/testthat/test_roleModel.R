#test_roleModel.R
library(testthat)
library(roleR)
context('several model runs work')

# works
test_that('model of 100 iters runs without crash with filled timeseries', {
  model <- roleR::dummyModel(R=FALSE,run=TRUE,fill_ts=F,niters=100)
  model$local$abundance_indv
})

# works
test_that('model of 100 iters runs without crash', {
  model <- roleR::dummyModel(R=FALSE,run=TRUE,fill_ts=FALSE,niters=100,print=TRUE)
})

# works
test_that('model of 1000 iters runs without crash', {
  model <- dummyModel(R=FALSE,run=TRUE,niters=1000)
})

# works
test_that('model of 1000 iters runs without crash when not doing speciation or immigration', {
  library(roleR)
  model <- roleR::dummyModel(R=FALSE,run=TRUE,niters=100,no_speciation=T,no_dispersal=T,print=F)
})