library(testthat)
library(roleR)

# works
test_that('model of 1000 iters is created successfully', {
  model <- roleR::dummyModel(R=F,run=T,fill_ts=F,niters=1000)
})

# works
test_that('model of 1000 iters is created successfully', {
  model <- roleR::dummyModel(R=F,run=F,fill_ts=F,niters=100)
  model <- roleR::iterSim(model,100,10,F)
})

# works
test_that('model of 1000 iters is created successfully', {
  model1 <- roleR::dummyModel(R=F,run=F,fill_ts=F,niters=10000)
  model1 <- roleR::iterSim(model1,10000,10,F)
  
  model2 <- roleR::dummyModel(R=F,run=T,fill_ts=F,niters=1000)
  model2 <- roleR::iterSim(model2,1000,10,F)
  
  model3 <- roleR::dummyModel(R=F,run=F,fill_ts=F,niters=10)
  model3$params$species_meta
  model3 <- roleR::iterSim(model3,1000,10,F)
  
})