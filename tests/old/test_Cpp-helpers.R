# test Cpp helpers

library(testthat)
library(roleR)

context('Cpp helpers work')

# works
test_that('copyData copys and copys deeply', {
  # script to run model on command line for testing
  library(roleR)
  model <- dummyModel(R=FALSE,run=FALSE)
  
  d0 <- model$copyData(1)

  iterSim(model,100,100,FALSE)
  
  d1 <- model$copyData(1)
  
  expect_false(sum(d0$local$abundance_sp == d1$local$abundance_sp) == 1000)
})